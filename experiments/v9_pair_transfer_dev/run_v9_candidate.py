from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
CLEAN_DIR = HERE / "clean_pair_poc"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(CLEAN_DIR))

from nabu_protein.v83 import NabuV83Model
from run_v9_pair_transfer_ired import (
    b2_predictions,
    derive_reference,
    evaluate,
    evaluate_by_mutation_count,
    fit_pair_factorization,
    load_ired,
    mutation_set,
    phase2_abstention_predictions,
    prediction_hash,
    score_v9,
)
from run_v9_1_degree_normalized import rank_hash, score_v91
from run_clean_pair_poc import (
    N_FOLDS,
    build_clean_pair_table,
    fit_node_additive,
    predict_node_additive,
)
from run_clean_pair_calibration import fit_affine


VERSION = "NABU_V9_CLEAN_CALIBRATED_PAIR_EXPANSION_DEV_V1"


def build_singleton_map(fit_frame: pd.DataFrame, reference: str):
    wt_target = None
    singleton_target = {}
    for sequence, target in zip(fit_frame["sequence"], fit_frame["target"]):
        mutations = mutation_set(sequence, reference)
        if len(mutations) == 0:
            if wt_target is not None:
                raise RuntimeError("Expected exactly one WT row.")
            wt_target = float(target)
        elif len(mutations) == 1:
            mutation = mutations[0]
            if mutation in singleton_target:
                raise RuntimeError(f"Duplicate singleton identity: {mutation}")
            singleton_target[mutation] = float(target)

    if wt_target is None:
        raise RuntimeError("Missing WT row.")

    clean_main = {
        mutation: float(target - wt_target)
        for mutation, target in singleton_target.items()
    }
    return float(wt_target), singleton_target, clean_main


def fit_global_oof_calibration(clean: pd.DataFrame):
    rows = []
    for fold in range(N_FOLDS):
        train = clean[clean["fold"] != fold]
        hold = clean[clean["fold"] == fold]

        model = fit_node_additive(
            list(zip(train["left"], train["right"])),
            train["clean_pair_residual"].to_numpy(dtype=float),
            permute=False,
            fold=fold,
        )
        train_nodes = set(model["effect"])

        for row in hold.itertuples(index=False):
            if row.left not in train_nodes or row.right not in train_nodes:
                continue
            raw = predict_node_additive((row.left, row.right), model)
            if raw is None:
                continue
            rows.append(
                {
                    "candidate_id": row.candidate_id,
                    "fold": int(fold),
                    "raw_prediction": float(raw),
                    "target": float(row.clean_pair_residual),
                }
            )

    oof = pd.DataFrame(rows)
    if len(oof) < 2:
        raise RuntimeError("Insufficient OOF rows for global calibration.")

    calibration = fit_affine(
        oof["raw_prediction"].to_numpy(dtype=float),
        oof["target"].to_numpy(dtype=float),
    )
    return calibration, oof


def exact_pair_map(clean: pd.DataFrame):
    result = {}
    for row in clean.itertuples(index=False):
        key = tuple(sorted((row.left, row.right)))
        if key in result:
            raise RuntimeError(f"Duplicate clean pair: {key}")
        result[key] = float(row.clean_pair_residual)
    return result


def calibrated_pair_prediction(pair, relation_model, calibration):
    raw = predict_node_additive(pair, relation_model)
    if raw is None:
        return None
    return float(calibration["alpha"] + calibration["beta"] * raw)


def score_v9_clean_candidate(
    mutations,
    *,
    wt_target,
    clean_main,
    exact_pairs,
    relation_model,
    calibration,
    b2_fallback,
):
    if not all(mutation in clean_main for mutation in mutations):
        return float(b2_fallback), {
            "complete_clean_main": False,
            "used_b2_fallback": True,
            "exact_clean_pair_count": 0,
            "calibrated_unseen_pair_count": 0,
            "unresolved_pair_count": (
                int(math.comb(len(mutations), 2))
                if len(mutations) >= 2
                else 0
            ),
            "pair_sum": 0.0,
        }

    score = float(wt_target + sum(clean_main[m] for m in mutations))
    exact_count = 0
    predicted_count = 0
    unresolved_count = 0
    pair_sum = 0.0

    for left, right in combinations(mutations, 2):
        key = tuple(sorted((left, right)))
        if key in exact_pairs:
            contribution = float(exact_pairs[key])
            exact_count += 1
        else:
            contribution = calibrated_pair_prediction(
                (left, right),
                relation_model,
                calibration,
            )
            if contribution is None:
                contribution = 0.0
                unresolved_count += 1
            else:
                predicted_count += 1

        pair_sum += float(contribution)

    score += pair_sum
    return float(score), {
        "complete_clean_main": True,
        "used_b2_fallback": False,
        "exact_clean_pair_count": int(exact_count),
        "calibrated_unseen_pair_count": int(predicted_count),
        "unresolved_pair_count": int(unresolved_count),
        "pair_sum": float(pair_sum),
    }


def summarize_scale(values):
    values = np.asarray(values, dtype=float)
    return {
        "mean": float(np.mean(values)),
        "std": float(np.std(values)),
        "min": float(np.min(values)),
        "max": float(np.max(values)),
    }


def run(source_gz: Path, output_dir: Path):
    source_sha, fit_frame, validation_frame, test_frame = load_ired(source_gz)
    reference = derive_reference(fit_frame["sequence"].tolist())

    fit_mutations = [
        mutation_set(sequence, reference)
        for sequence in fit_frame["sequence"]
    ]
    test_mutations = [
        mutation_set(sequence, reference)
        for sequence in test_frame["sequence"]
    ]
    fit_target = fit_frame["target"].to_numpy(dtype=float)
    test_target = test_frame["target"].to_numpy(dtype=float)
    test_ids = test_frame["sequence"].tolist()

    v83_model = NabuV83Model().fit(
        mutation_sets=fit_mutations,
        labels=fit_target,
        candidate_ids=fit_frame["sequence"].tolist(),
    )

    phase2_prediction, strict_scoreable = phase2_abstention_predictions(
        v83_model,
        test_mutations,
        test_ids,
    )
    b2_prediction, b2_all_main = b2_predictions(
        v83_model,
        test_mutations,
    )

    old_factor = fit_pair_factorization(
        v83_model,
        fit_mutations,
        fit_target,
        permute_targets=False,
    )
    v90_prediction = np.asarray(
        [
            score_v9(mutations, v83_model, old_factor)[0]
            for mutations in test_mutations
        ],
        dtype=float,
    )
    v91_prediction = np.asarray(
        [
            score_v91(mutations, v83_model, old_factor)[0]
            for mutations in test_mutations
        ],
        dtype=float,
    )

    clean_metadata, clean = build_clean_pair_table(source_gz)
    wt_target, singleton_target, clean_main = build_singleton_map(
        fit_frame,
        reference,
    )
    clean_exact_pairs = exact_pair_map(clean)
    calibration, calibration_oof = fit_global_oof_calibration(clean)

    relation_model = fit_node_additive(
        list(zip(clean["left"], clean["right"])),
        clean["clean_pair_residual"].to_numpy(dtype=float),
        permute=False,
        fold=0,
    )

    rows = [
        score_v9_clean_candidate(
            mutations,
            wt_target=wt_target,
            clean_main=clean_main,
            exact_pairs=clean_exact_pairs,
            relation_model=relation_model,
            calibration=calibration,
            b2_fallback=b2_value,
        )
        for mutations, b2_value in zip(test_mutations, b2_prediction)
    ]
    v9_prediction = np.asarray([value for value, _ in rows], dtype=float)
    diag = pd.DataFrame([item for _, item in rows])

    replay_rows = [
        score_v9_clean_candidate(
            mutations,
            wt_target=wt_target,
            clean_main=clean_main,
            exact_pairs=clean_exact_pairs,
            relation_model=relation_model,
            calibration=calibration,
            b2_fallback=b2_value,
        )
        for mutations, b2_value in zip(test_mutations, b2_prediction)
    ]
    replay_prediction = np.asarray(
        [value for value, _ in replay_rows],
        dtype=float,
    )

    arms = {
        "PHASE2_V8_3_ABSTENTION": phase2_prediction,
        "B2_MAIN_ONLY": b2_prediction,
        "V9_0_RAW_PAIR_TRANSFER": v90_prediction,
        "V9_1_DEGREE_NORMALIZED_TRANSFER": v91_prediction,
        "V9_CLEAN_CALIBRATED_PAIR_EXPANSION": v9_prediction,
    }

    if not all(np.isfinite(values).all() for values in arms.values()):
        raise RuntimeError("Non-finite prediction detected.")

    metrics = {
        name: evaluate(test_target, values)
        for name, values in arms.items()
    }
    by_mutation_count = {
        name: evaluate_by_mutation_count(
            test_target,
            values,
            test_mutations,
        )
        for name, values in arms.items()
    }

    clean_mask = diag["complete_clean_main"].to_numpy(dtype=bool)
    if int(clean_mask.sum()) < 2:
        raise RuntimeError("Insufficient complete-clean-main test rows.")

    clean_region_metrics = {
        name: evaluate(test_target[clean_mask], values[clean_mask])
        for name, values in {
            "B2_MAIN_ONLY": b2_prediction,
            "V9_1_DEGREE_NORMALIZED_TRANSFER": v91_prediction,
            "V9_CLEAN_CALIBRATED_PAIR_EXPANSION": v9_prediction,
        }.items()
    }

    target_std = float(np.std(test_target))
    v9_std = float(np.std(v9_prediction))
    v9_metrics = metrics["V9_CLEAN_CALIBRATED_PAIR_EXPANSION"]
    b2_metrics = metrics["B2_MAIN_ONLY"]
    v91_metrics = metrics["V9_1_DEGREE_NORMALIZED_TRANSFER"]

    count_metrics = by_mutation_count[
        "V9_CLEAN_CALIBRATED_PAIR_EXPANSION"
    ]
    positive_3_4_5 = all(
        str(count) in count_metrics
        and count_metrics[str(count)]["spearman"] is not None
        and count_metrics[str(count)]["spearman"] > 0.0
        for count in (3, 4, 5)
    )

    checks = {
        "full_spearman_beats_b2": bool(
            v9_metrics["spearman"] is not None
            and b2_metrics["spearman"] is not None
            and v9_metrics["spearman"] > b2_metrics["spearman"]
        ),
        "full_spearman_beats_v9_1": bool(
            v9_metrics["spearman"] is not None
            and v91_metrics["spearman"] is not None
            and v9_metrics["spearman"] > v91_metrics["spearman"]
        ),
        "top1_hits_at_least_b2": bool(
            v9_metrics["top1_percent_hits"]
            >= b2_metrics["top1_percent_hits"]
        ),
        "regret_not_worse_than_b2": bool(
            v9_metrics["normalized_regret_top1pct"]
            <= b2_metrics["normalized_regret_top1pct"]
        ),
        "clean_region_spearman_beats_b2": bool(
            clean_region_metrics["V9_CLEAN_CALIBRATED_PAIR_EXPANSION"][
                "spearman"
            ]
            > clean_region_metrics["B2_MAIN_ONLY"]["spearman"]
        ),
        "positive_spearman_counts_3_4_5": bool(positive_3_4_5),
        "prediction_scale_within_3x_target_std": bool(
            np.isfinite(v9_std)
            and target_std > 0.0
            and v9_std <= 3.0 * target_std
        ),
        "deterministic_rank_replay": bool(
            rank_hash(test_ids, v9_prediction)
            == rank_hash(test_ids, replay_prediction)
        ),
    }

    predictions = pd.DataFrame(
        {
            "candidate_id": test_ids,
            "mutation_count": [len(x) for x in test_mutations],
            "target": test_target,
            "strict_v83_scoreable": strict_scoreable,
            "b2_all_main_supported": b2_all_main,
            **arms,
            "v9_complete_clean_main": diag["complete_clean_main"],
            "v9_used_b2_fallback": diag["used_b2_fallback"],
            "v9_exact_clean_pair_count": diag["exact_clean_pair_count"],
            "v9_calibrated_unseen_pair_count": diag[
                "calibrated_unseen_pair_count"
            ],
            "v9_unresolved_pair_count": diag["unresolved_pair_count"],
            "v9_pair_sum": diag["pair_sum"],
        }
    )

    result = {
        "version": VERSION,
        "status": "POST_PHASE2_DEVELOPMENT_DIAGNOSTIC",
        "phase3_opened": False,
        "scientific_claim": False,
        "source_sha256": source_sha,
        "frozen_v83_core": (
            "c8afdcd7698231d95a39ef13a3fe22b6e3f507c5"
        ),
        "split_counts": {
            "fit": int(len(fit_frame)),
            "validation_unused_for_fit": int(len(validation_frame)),
            "test_development": int(len(test_frame)),
        },
        "reference_sha256": hashlib.sha256(
            reference.encode("utf-8")
        ).hexdigest(),
        "reference_length": int(len(reference)),
        "clean_evidence": {
            "wt_target": float(wt_target),
            "singleton_count": int(len(singleton_target)),
            "clean_main_count": int(len(clean_main)),
            "clean_pair_count": int(len(clean)),
            "exact_clean_pair_count": int(len(clean_exact_pairs)),
            "relation_node_count": int(len(relation_model["effect"])),
            "calibration_oof_count": int(len(calibration_oof)),
            "calibration_alpha": float(calibration["alpha"]),
            "calibration_beta": float(calibration["beta"]),
            "calibration_rank": int(calibration["rank"]),
        },
        "coverage": {
            "strict_v83_scoreable": int(strict_scoreable.sum()),
            "b2_all_main_supported": int(b2_all_main.sum()),
            "complete_clean_main": int(clean_mask.sum()),
            "b2_fallback_rows": int(diag["used_b2_fallback"].sum()),
            "rows_with_exact_clean_pairs": int(
                (diag["exact_clean_pair_count"] > 0).sum()
            ),
            "rows_with_calibrated_unseen_pairs": int(
                (diag["calibrated_unseen_pair_count"] > 0).sum()
            ),
            "exact_clean_pair_contributions": int(
                diag["exact_clean_pair_count"].sum()
            ),
            "calibrated_unseen_pair_contributions": int(
                diag["calibrated_unseen_pair_count"].sum()
            ),
            "unresolved_pair_contributions": int(
                diag["unresolved_pair_count"].sum()
            ),
        },
        "metrics": metrics,
        "complete_clean_main_region_metrics": clean_region_metrics,
        "by_mutation_count": by_mutation_count,
        "prediction_scale": {
            name: summarize_scale(values)
            for name, values in arms.items()
        },
        "target_scale": summarize_scale(test_target),
        "prediction_sha256": {
            name: prediction_hash(test_ids, values)
            for name, values in arms.items()
        },
        "rank_sha256": {
            name: rank_hash(test_ids, values)
            for name, values in arms.items()
        },
        "replay_rank_sha256": rank_hash(
            test_ids,
            replay_prediction,
        ),
        "decision_checks": checks,
        "development_decision": (
            "V9_CLEAN_CALIBRATED_PAIR_EXPANSION_RETAIN"
            if all(checks.values())
            else "V9_CLEAN_CALIBRATED_PAIR_EXPANSION_REJECT"
        ),
        "interpretation_boundary": (
            "IRED test targets were already revealed during Phase 2D. "
            "This is development evidence only. A retained V9 candidate "
            "still requires a new untouched external benchmark."
        ),
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    predictions.to_csv(
        output_dir / "V9_CANDIDATE_PREDICTIONS.csv",
        index=False,
    )
    (output_dir / "V9_CANDIDATE_RESULTS.json").write_text(
        json.dumps(result, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2))
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source_gz")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    run(Path(args.source_gz), Path(args.out))


if __name__ == "__main__":
    main()
