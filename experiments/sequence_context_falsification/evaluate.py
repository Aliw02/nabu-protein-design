from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

SEED = 161
EPS = 1e-12
PERMUTATIONS = 1000
LOW_FEATURES = [f"low_feature_{i:02d}" for i in range(15)]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def finite_spearman(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) < 2 or np.std(a) <= EPS or np.std(b) <= EPS:
        return None
    value = float(spearmanr(a, b).statistic)
    return value if np.isfinite(value) else None


def finite_pearson(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) < 2 or np.std(a) <= EPS or np.std(b) <= EPS:
        return None
    value = float(pearsonr(a, b).statistic)
    return value if np.isfinite(value) else None


def residualize(X, y):
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    design = np.column_stack([np.ones(len(X), dtype=float), X])
    beta, _, rank, singular = np.linalg.lstsq(design, y, rcond=None)
    fitted = design @ beta
    residual = y - fitted
    return residual, {
        "design_rank": int(rank),
        "design_columns": int(design.shape[1]),
        "residual_std": float(np.std(residual)),
        "fitted_std": float(np.std(fitted)),
        "singular_value_min": float(np.min(singular)),
        "singular_value_max": float(np.max(singular)),
    }


def affine_eval(target, score):
    target = np.asarray(target, dtype=float)
    score = np.asarray(score, dtype=float)
    design = np.column_stack(
        [np.ones(len(score), dtype=float), score]
    )
    beta, *_ = np.linalg.lstsq(design, target, rcond=None)
    pred = design @ beta
    return {
        "intercept": float(beta[0]),
        "slope": float(beta[1]),
        "rmse": float(np.sqrt(np.mean((target - pred) ** 2))),
        "target_std": float(np.std(target)),
        "normalized_rmse": (
            None
            if np.std(target) <= EPS
            else float(
                np.sqrt(np.mean((target - pred) ** 2))
                / np.std(target)
            )
        ),
        "evaluation_only": True,
    }


def evaluate_cell(cell, dataset, order):
    X = cell[LOW_FEATURES].to_numpy(dtype=float)
    target = cell["residual_ho"].to_numpy(dtype=float)
    context = cell["esm_context_mean"].to_numpy(dtype=float)
    single = cell["esm_single_mean"].to_numpy(dtype=float)

    if len(cell) != 192:
        raise RuntimeError(
            f"{dataset} order {order}: expected 192 rows, got {len(cell)}"
        )

    target_resid, target_diag = residualize(X, target)
    context_resid, context_diag = residualize(X, context)
    single_resid, single_diag = residualize(X, single)

    partial_context = finite_spearman(target_resid, context_resid)
    partial_single = finite_spearman(target_resid, single_resid)

    rng = np.random.default_rng(
        SEED + (0 if dataset == "GB1" else 10000) + int(order)
    )
    null_values = []
    for _ in range(PERMUTATIONS):
        permuted = rng.permutation(context)
        perm_resid, _ = residualize(X, permuted)
        null_values.append(
            finite_spearman(target_resid, perm_resid)
        )
    null = np.asarray(null_values, dtype=float)
    if not np.isfinite(null).all():
        raise RuntimeError("Permutation null contains non-finite values.")

    observed_abs = abs(float(partial_context))
    empirical_p = float(
        (1 + np.sum(np.abs(null) >= observed_abs))
        / (PERMUTATIONS + 1)
    )

    metrics = {
        "dataset": dataset,
        "order": int(order),
        "count": int(len(cell)),
        "raw_context": {
            "spearman": finite_spearman(target, context),
            "pearson": finite_pearson(target, context),
            "affine_fit": affine_eval(target, context),
        },
        "raw_single_control": {
            "spearman": finite_spearman(target, single),
            "pearson": finite_pearson(target, single),
        },
        "low_order_controlled": {
            "partial_context_spearman": partial_context,
            "partial_single_spearman": partial_single,
            "abs_context_minus_abs_single": float(
                abs(float(partial_context))
                - abs(float(partial_single))
            ),
            "target_residualization": target_diag,
            "context_residualization": context_diag,
            "single_residualization": single_diag,
        },
        "permutation": {
            "count": PERMUTATIONS,
            "empirical_two_sided_p": empirical_p,
            "null_mean": float(np.mean(null)),
            "null_std": float(np.std(null)),
            "null_abs_p95": float(
                np.quantile(np.abs(null), 0.95)
            ),
            "null_abs_p99": float(
                np.quantile(np.abs(null), 0.99)
            ),
            "observed_abs_partial_context_spearman": observed_abs,
        },
        "score_distribution": {
            "context_mean": float(np.mean(context)),
            "context_std": float(np.std(context)),
            "context_min": float(np.min(context)),
            "context_max": float(np.max(context)),
            "single_mean": float(np.mean(single)),
            "single_std": float(np.std(single)),
            "target_residual_mean": float(np.mean(target)),
            "target_residual_std": float(np.std(target)),
        },
    }
    return metrics, null


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--visible", required=True)
    parser.add_argument("--scores", required=True)
    parser.add_argument("--reveal", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    visible_dir = Path(args.visible)
    scores_dir = Path(args.scores)
    reveal_dir = Path(args.reveal)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    selected = pd.read_csv(
        visible_dir / "SELECTED_CANDIDATES.csv"
    )
    scores = pd.read_csv(
        scores_dir / "ESM_CANDIDATE_SCORES.csv"
    )
    reveal = pd.read_csv(
        reveal_dir / "REVEAL_TARGETS.csv"
    )

    keys = ["dataset", "order", "candidate_id"]
    for frame, name in (
        (selected, "selected"),
        (scores, "scores"),
        (reveal, "reveal"),
    ):
        if frame.duplicated(keys).any():
            raise RuntimeError(f"Duplicate keys in {name}.")

    joined = (
        selected.merge(scores, on=keys, how="inner", validate="one_to_one")
        .merge(reveal, on=keys, how="inner", validate="one_to_one")
    )
    if len(joined) != 768:
        raise RuntimeError(
            f"Joined evidence row mismatch: {len(joined)}"
        )

    model_manifest = json.loads(
        (scores_dir / "MODEL_MANIFEST.json").read_text()
    )
    scoring_diag = json.loads(
        (scores_dir / "SCORING_DIAGNOSTICS.json").read_text()
    )
    source_manifest = json.loads(
        (visible_dir / "SOURCE_MANIFEST.json").read_text()
    )

    required_checks = {
        "model_revision_exact": (
            model_manifest["resolved_revision"]
            == "c731040fcd8d73dceaa04b0a8e6329b345b0f5df"
        ),
        "model_weights_hash_exact": (
            model_manifest["files"]["model.safetensors"]["sha256"]
            == "24c5fa474c48f3b754b86efe752d5f189d2bcd88190fa2270fc92b2ef3034189"
        ),
        "deterministic_replay_exact": bool(
            scoring_diag["determinism"]["exact"]
            and scoring_diag["determinism"][
                "max_abs_probability_score_difference"
            ] == 0.0
        ),
        "nucb_not_consumed": (
            source_manifest["nucb_consumed"] is False
        ),
        "phase3_closed": (
            source_manifest["phase3_opened"] is False
        ),
        "selected_row_count_exact": len(selected) == 768,
        "score_row_count_exact": len(scores) == 768,
        "reveal_row_count_exact": len(reveal) == 768,
    }
    all_checks = bool(all(required_checks.values()))

    cell_metrics = {}
    permutation_rows = []

    for dataset in ("TrpB", "GB1"):
        cell_metrics[dataset] = {}
        for order in (3, 4):
            cell = joined[
                (joined["dataset"] == dataset)
                & (joined["order"] == order)
            ].copy()
            metrics, null = evaluate_cell(
                cell,
                dataset,
                order,
            )
            cell_metrics[dataset][str(order)] = metrics
            for i, value in enumerate(null):
                permutation_rows.append(
                    {
                        "dataset": dataset,
                        "order": order,
                        "permutation_index": i,
                        "partial_spearman": float(value),
                    }
                )

    quad = {
        dataset: cell_metrics[dataset]["4"]
        for dataset in ("TrpB", "GB1")
    }

    quad_conditions = {}
    for dataset, block in quad.items():
        partial = float(
            block["low_order_controlled"][
                "partial_context_spearman"
            ]
        )
        single = float(
            block["low_order_controlled"][
                "partial_single_spearman"
            ]
        )
        p = float(
            block["permutation"]["empirical_two_sided_p"]
        )
        quad_conditions[dataset] = {
            "abs_partial_at_least_0_15": abs(partial) >= 0.15,
            "permutation_p_at_most_0_01": p <= 0.01,
            "beats_single_abs_by_0_05": (
                abs(partial) >= abs(single) + 0.05
            ),
        }

    trpb_partial = float(
        quad["TrpB"]["low_order_controlled"][
            "partial_context_spearman"
        ]
    )
    gb1_partial = float(
        quad["GB1"]["low_order_controlled"][
            "partial_context_spearman"
        ]
    )
    same_sign = bool(
        np.sign(trpb_partial) == np.sign(gb1_partial)
        and np.sign(trpb_partial) != 0
    )

    supported = bool(
        all_checks
        and same_sign
        and all(
            all(conditions.values())
            for conditions in quad_conditions.values()
        )
    )

    weak = {}
    for dataset, block in quad.items():
        partial = float(
            block["low_order_controlled"][
                "partial_context_spearman"
            ]
        )
        p = float(
            block["permutation"]["empirical_two_sided_p"]
        )
        weak[dataset] = bool(
            abs(partial) < 0.08 or p > 0.10
        )

    if supported:
        decision = (
            "PRETRAINED_SEQUENCE_CONTEXT_SIGNAL_SUPPORTED_FOR_QUADS"
        )
        next_action = (
            "STOP_DIAGNOSTICS_PREREGISTER_ONE_SEQUENCE_CONTEXT_ARCHITECTURE"
        )
    elif all(weak.values()):
        decision = "PRETRAINED_SEQUENCE_CONTEXT_SCALAR_REJECTED"
        next_action = (
            "DO_NOT_BUILD_FROM_THIS_SCALAR_SCORE"
        )
    else:
        decision = "PRETRAINED_SEQUENCE_CONTEXT_SCALAR_UNRESOLVED"
        next_action = "NO_ARCHITECTURE_PRESERVE_EVIDENCE"

    decision_block = {
        "version": "NABU_SEQUENCE_CONTEXT_FALSIFICATION_V1",
        "decision": decision,
        "next_action": next_action,
        "root_cause_status": "ROOT_CAUSE_NOT_YET_IDENTIFIED",
        "quad_gate": quad_conditions,
        "quad_same_sign": same_sign,
        "quad_partial_context_spearman": {
            "TrpB": trpb_partial,
            "GB1": gb1_partial,
        },
        "quad_weak_rejection_condition": weak,
        "required_checks": required_checks,
        "all_required_checks_pass": all_checks,
        "nucb_consumed": False,
        "phase3_opened": False,
    }

    joined.to_csv(
        out / "JOINED_EVIDENCE.csv",
        index=False,
    )
    pd.DataFrame(permutation_rows).to_csv(
        out / "PERMUTATION_NULL.csv",
        index=False,
    )
    (out / "CELL_METRICS.json").write_text(
        json.dumps(cell_metrics, indent=2),
        encoding="utf-8",
    )
    (out / "DECISION.json").write_text(
        json.dumps(decision_block, indent=2),
        encoding="utf-8",
    )

    evaluation_diag = {
        "joined_rows": int(len(joined)),
        "cells": {
            f"{dataset}_order_{order}": int(
                np.sum(
                    (joined["dataset"] == dataset)
                    & (joined["order"] == order)
                )
            )
            for dataset in ("TrpB", "GB1")
            for order in (3, 4)
        },
        "low_feature_columns": LOW_FEATURES,
        "permutations_per_cell": PERMUTATIONS,
        "seed": SEED,
        "source_manifest_sha256": sha256_file(
            visible_dir / "SOURCE_MANIFEST.json"
        ),
        "selected_candidates_sha256": sha256_file(
            visible_dir / "SELECTED_CANDIDATES.csv"
        ),
        "esm_candidate_scores_sha256": sha256_file(
            scores_dir / "ESM_CANDIDATE_SCORES.csv"
        ),
        "reveal_targets_sha256": sha256_file(
            reveal_dir / "REVEAL_TARGETS.csv"
        ),
    }
    (out / "EVALUATION_DIAGNOSTICS.json").write_text(
        json.dumps(evaluation_diag, indent=2),
        encoding="utf-8",
    )

    hashes = {
        p.name: sha256_file(p)
        for p in sorted(out.iterdir())
        if p.is_file() and p.name != "OUTPUT_HASHES.json"
    }
    (out / "OUTPUT_HASHES.json").write_text(
        json.dumps(hashes, indent=2),
        encoding="utf-8",
    )

    print(json.dumps(decision_block, indent=2))


if __name__ == "__main__":
    main()
