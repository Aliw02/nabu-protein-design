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
EXPECTED_MODEL_REVISION = "c731040fcd8d73dceaa04b0a8e6329b345b0f5df"
EXPECTED_MODEL_SHA256 = (
    "24c5fa474c48f3b754b86efe752d5f189d2bcd88190fa2270fc92b2ef3034189"
)
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
    design = np.column_stack(
        [np.ones(len(X), dtype=float), X]
    )
    beta, _, rank, singular = np.linalg.lstsq(
        design,
        y,
        rcond=None,
    )
    fitted = design @ beta
    return (
        y - fitted,
        {
            "rank": int(rank),
            "column_count": int(design.shape[1]),
            "min_singular_value": float(np.min(singular)),
            "max_singular_value": float(np.max(singular)),
            "beta": [float(x) for x in beta],
        },
    )


def permutation_test(X, y_residual, score, observed):
    rng = np.random.default_rng(SEED)
    null = np.empty(PERMUTATIONS, dtype=float)

    for i in range(PERMUTATIONS):
        shuffled = rng.permutation(score)
        score_resid, _ = residualize(X, shuffled)
        rho = finite_spearman(y_residual, score_resid)
        if rho is None:
            raise RuntimeError("Permutation produced undefined Spearman.")
        null[i] = rho

    p = float(
        (1 + np.sum(np.abs(null) >= abs(float(observed))))
        / (PERMUTATIONS + 1)
    )

    return null, {
        "permutations": PERMUTATIONS,
        "empirical_two_sided_p": p,
        "null_mean": float(np.mean(null)),
        "null_std": float(np.std(null)),
        "null_abs_rho_p95": float(np.quantile(np.abs(null), 0.95)),
        "null_abs_rho_p99": float(np.quantile(np.abs(null), 0.99)),
    }


def verify_hash_manifest(directory: Path, manifest_name: str):
    manifest_path = directory / manifest_name
    if not manifest_path.exists():
        raise RuntimeError(f"Missing hash manifest: {manifest_path}")
    manifest = json.loads(manifest_path.read_text())
    checks = {}
    for name, expected in manifest.items():
        path = directory / name
        if not path.exists():
            raise RuntimeError(f"Hashed file missing: {path}")
        observed = sha256_file(path)
        checks[name] = {
            "expected": expected,
            "observed": observed,
            "match": bool(observed == expected),
        }
        if observed != expected:
            raise RuntimeError(f"Hash mismatch for {path}")
    return checks


def evaluate_cell(cell, permutation_rows):
    X = cell[LOW_FEATURES].to_numpy(dtype=float)
    y = cell["residual_ho"].to_numpy(dtype=float)
    mobius = cell["model_ho"].to_numpy(dtype=float)
    context = cell["context_mean_control"].to_numpy(dtype=float)
    single = cell["single_mean_control"].to_numpy(dtype=float)

    y_resid, y_diag = residualize(X, y)
    mobius_resid, mobius_diag = residualize(X, mobius)
    context_resid, context_diag = residualize(X, context)
    single_resid, single_diag = residualize(X, single)

    partial_mobius = finite_spearman(y_resid, mobius_resid)
    partial_context = finite_spearman(y_resid, context_resid)
    partial_single = finite_spearman(y_resid, single_resid)

    if partial_mobius is None:
        raise RuntimeError("Undefined partial Mobius Spearman.")

    null, null_summary = permutation_test(
        X,
        y_resid,
        mobius,
        partial_mobius,
    )

    dataset = str(cell["dataset"].iloc[0])
    order = int(cell["order"].iloc[0])
    for i, rho in enumerate(null):
        permutation_rows.append(
            {
                "dataset": dataset,
                "order": order,
                "permutation_index": i,
                "partial_mobius_spearman": float(rho),
            }
        )

    return {
        "count": int(len(cell)),
        "selection_rank_min": int(cell["selection_rank"].min()),
        "selection_rank_max": int(cell["selection_rank"].max()),
        "raw": {
            "mobius_spearman": finite_spearman(y, mobius),
            "mobius_pearson": finite_pearson(y, mobius),
            "context_mean_spearman": finite_spearman(y, context),
            "context_mean_pearson": finite_pearson(y, context),
            "single_mean_spearman": finite_spearman(y, single),
            "single_mean_pearson": finite_pearson(y, single),
            "bio_ho_std": float(np.std(y)),
            "model_ho_std": float(np.std(mobius)),
        },
        "controlled": {
            "partial_mobius_spearman": partial_mobius,
            "partial_context_spearman": partial_context,
            "partial_single_spearman": partial_single,
            "abs_mobius_minus_abs_context": float(
                abs(partial_mobius) - abs(partial_context)
            ),
            "abs_mobius_minus_abs_single": float(
                abs(partial_mobius) - abs(partial_single)
            ),
        },
        "permutation": null_summary,
        "residualization": {
            "bio_ho": y_diag,
            "model_ho": mobius_diag,
            "context_mean": context_diag,
            "single_mean": single_diag,
        },
    }


def quad_gate(cell):
    c = cell["controlled"]
    p = cell["permutation"]["empirical_two_sided_p"]
    return {
        "abs_partial_mobius_at_least_0_15": bool(
            abs(c["partial_mobius_spearman"]) >= 0.15
        ),
        "permutation_p_at_most_0_01": bool(p <= 0.01),
        "beats_context_by_abs_0_05": bool(
            abs(c["partial_mobius_spearman"])
            >= abs(c["partial_context_spearman"]) + 0.05
        ),
        "beats_single_by_abs_0_05": bool(
            abs(c["partial_mobius_spearman"])
            >= abs(c["partial_single_spearman"]) + 0.05
        ),
    }


def weak_rejection(cell):
    rho = abs(cell["controlled"]["partial_mobius_spearman"])
    p = cell["permutation"]["empirical_two_sided_p"]
    return bool(rho < 0.08 or p > 0.10)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--visible", required=True)
    parser.add_argument("--scores", required=True)
    parser.add_argument("--reveal", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    visible_dir = Path(args.visible)
    score_dir = Path(args.scores)
    reveal_dir = Path(args.reveal)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    visible_hash_checks = verify_hash_manifest(
        visible_dir,
        "VISIBLE_HASHES.json",
    )
    score_hash_checks = verify_hash_manifest(
        score_dir,
        "SCORE_HASHES.json",
    )
    reveal_hash_checks = verify_hash_manifest(
        reveal_dir,
        "REVEAL_HASHES.json",
    )

    selected = pd.read_csv(
        visible_dir / "SELECTED_CANDIDATES.csv"
    )
    scores = pd.read_csv(
        score_dir / "CANDIDATE_MODEL_SCORES.csv"
    )
    reveal = pd.read_csv(
        reveal_dir / "REVEAL_TARGETS.csv"
    )

    expected = 2 * 2 * 192
    for label, frame in (
        ("selected", selected),
        ("scores", scores),
        ("reveal", reveal),
    ):
        if len(frame) != expected:
            raise RuntimeError(
                f"{label} row count mismatch: {len(frame)} != {expected}"
            )

    keys = ["dataset", "order", "candidate_id"]
    if selected.duplicated(keys).any():
        raise RuntimeError("Selected candidate keys are not unique.")
    if scores.duplicated(keys).any():
        raise RuntimeError("Score candidate keys are not unique.")
    if reveal.duplicated(keys).any():
        raise RuntimeError("Reveal candidate keys are not unique.")

    score_forbidden = {
        "raw_target",
        "transformed_target",
        "exact_o2",
        "residual_ho",
    }
    if score_forbidden & set(scores.columns):
        raise RuntimeError(
            "Scoring artifact contains forbidden reveal columns."
        )

    joined = selected.merge(
        scores,
        on=["dataset", "order", "candidate_id", "selection_rank", "selection_hash"],
        how="inner",
        validate="one_to_one",
    ).merge(
        reveal,
        on=keys,
        how="inner",
        validate="one_to_one",
    )
    if len(joined) != expected:
        raise RuntimeError("Joined evidence count mismatch.")

    for feature in LOW_FEATURES:
        if feature not in joined.columns:
            raise RuntimeError(f"Missing low-order feature: {feature}")

    if int(joined["selection_rank"].min()) != 192:
        raise RuntimeError("Selection rank start is not 192.")
    if int(joined["selection_rank"].max()) != 383:
        raise RuntimeError("Selection rank end is not 383.")

    source_manifest = json.loads(
        (visible_dir / "SOURCE_MANIFEST.json").read_text()
    )
    model_manifest = json.loads(
        (score_dir / "MODEL_MANIFEST.json").read_text()
    )
    runtime_manifest = json.loads(
        (score_dir / "RUNTIME_MANIFEST.json").read_text()
    )
    scoring_diag = json.loads(
        (score_dir / "SCORING_DIAGNOSTICS.json").read_text()
    )

    model_files = model_manifest.get("files", {})
    weights = model_files.get("model.safetensors", {})

    integrity = {
        "model_revision_exact": bool(
            model_manifest.get("resolved_revision")
            == EXPECTED_MODEL_REVISION
        ),
        "model_weights_hash_exact": bool(
            weights.get("sha256") == EXPECTED_MODEL_SHA256
        ),
        "deterministic_replay_exact": bool(
            scoring_diag.get("determinism", {}).get("exact") is True
        ),
        "fresh_holdout_ranks_exact": bool(
            int(joined["selection_rank"].min()) == 192
            and int(joined["selection_rank"].max()) == 383
        ),
        "nucb_not_consumed": bool(
            source_manifest.get("nucb_consumed") is False
            and runtime_manifest.get("nucb_consumed") is False
            and scoring_diag.get("nucb_consumed") is False
        ),
        "phase3_closed": bool(
            source_manifest.get("phase3_opened") is False
            and runtime_manifest.get("phase3_opened") is False
            and scoring_diag.get("phase3_opened") is False
        ),
        "visible_hashes_exact": bool(
            all(x["match"] for x in visible_hash_checks.values())
        ),
        "score_hashes_exact": bool(
            all(x["match"] for x in score_hash_checks.values())
        ),
        "reveal_hashes_exact": bool(
            all(x["match"] for x in reveal_hash_checks.values())
        ),
    }
    integrity["all_required_pass"] = bool(all(integrity.values()))

    permutation_rows = []
    cells = {}
    for (dataset, order), cell in joined.groupby(
        ["dataset", "order"],
        sort=True,
    ):
        key = f"{dataset}|{int(order)}"
        cells[key] = evaluate_cell(
            cell.reset_index(drop=True),
            permutation_rows,
        )

    trpb4 = cells["TrpB|4"]
    gb14 = cells["GB1|4"]

    tr_gate = quad_gate(trpb4)
    gb_gate = quad_gate(gb14)

    tr_rho = trpb4["controlled"]["partial_mobius_spearman"]
    gb_rho = gb14["controlled"]["partial_mobius_spearman"]
    same_sign = bool(
        tr_rho != 0.0
        and gb_rho != 0.0
        and np.sign(tr_rho) == np.sign(gb_rho)
    )

    tr_gate["same_sign_across_quad_landscapes"] = same_sign
    gb_gate["same_sign_across_quad_landscapes"] = same_sign
    tr_gate["integrity_all_required_pass"] = integrity["all_required_pass"]
    gb_gate["integrity_all_required_pass"] = integrity["all_required_pass"]

    tr_pass = bool(all(tr_gate.values()))
    gb_pass = bool(all(gb_gate.values()))
    supported = bool(tr_pass and gb_pass)

    tr_reject = weak_rejection(trpb4)
    gb_reject = weak_rejection(gb14)

    if supported:
        decision = "PRETRAINED_MOBIUS_SIGNAL_SUPPORTED_FOR_QUADS"
        next_action = "PREREGISTER_ONE_ARCHITECTURE_NO_NUCB_YET"
    elif tr_reject and gb_reject:
        decision = "PRETRAINED_MOBIUS_SCALAR_REJECTED"
        next_action = (
            "DO_NOT_BUILD_FROM_ESM_LIKELIHOOD_ALGEBRA;"
            "RICHER_REPRESENTATION_REQUIRES_NEW_HYPOTHESIS"
        )
    else:
        decision = "PRETRAINED_MOBIUS_SCALAR_UNRESOLVED"
        next_action = "DO_NOT_BUILD_ARCHITECTURE"

    metrics = {
        "version": "NABU_PRETRAINED_MOBIUS_FALSIFICATION_METRICS_V1",
        "cells": cells,
        "quad_gate": {
            "TrpB": tr_gate,
            "GB1": gb_gate,
            "TrpB_all_required": tr_pass,
            "GB1_all_required": gb_pass,
        },
    }
    (out / "CELL_METRICS.json").write_text(
        json.dumps(metrics, indent=2),
        encoding="utf-8",
    )

    decision_doc = {
        "version": "NABU_PRETRAINED_MOBIUS_FALSIFICATION_DECISION_V1",
        "decision": decision,
        "next_action": next_action,
        "root_cause_status": "ROOT_CAUSE_NOT_YET_IDENTIFIED",
        "TrpB_quad_partial_mobius_spearman": tr_rho,
        "GB1_quad_partial_mobius_spearman": gb_rho,
        "same_sign": same_sign,
        "TrpB_weak_rejection_condition": tr_reject,
        "GB1_weak_rejection_condition": gb_reject,
        "integrity_all_required_pass": integrity["all_required_pass"],
        "nucb_consumed": False,
        "phase3_opened": False,
    }
    (out / "DECISION.json").write_text(
        json.dumps(decision_doc, indent=2),
        encoding="utf-8",
    )

    diagnostics = {
        "integrity": integrity,
        "visible_hash_checks": visible_hash_checks,
        "score_hash_checks": score_hash_checks,
        "reveal_hash_checks": reveal_hash_checks,
        "joined_count": int(len(joined)),
        "cell_counts": {
            f"{dataset}|{int(order)}": int(len(cell))
            for (dataset, order), cell in joined.groupby(
                ["dataset", "order"]
            )
        },
        "selection_rank_min": int(joined["selection_rank"].min()),
        "selection_rank_max": int(joined["selection_rank"].max()),
        "permutation_count_per_cell": PERMUTATIONS,
    }
    (out / "EVALUATION_DIAGNOSTICS.json").write_text(
        json.dumps(diagnostics, indent=2),
        encoding="utf-8",
    )

    joined.to_csv(out / "JOINED_EVIDENCE.csv", index=False)
    pd.DataFrame(permutation_rows).to_csv(
        out / "PERMUTATION_NULL.csv",
        index=False,
    )

    hashes = {
        path.name: sha256_file(path)
        for path in sorted(out.iterdir())
        if path.is_file() and path.name != "OUTPUT_HASHES.json"
    }
    (out / "OUTPUT_HASHES.json").write_text(
        json.dumps(hashes, indent=2),
        encoding="utf-8",
    )

    print(json.dumps(decision_doc, indent=2))


if __name__ == "__main__":
    main()
