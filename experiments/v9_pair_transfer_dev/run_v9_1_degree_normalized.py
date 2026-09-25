from __future__ import annotations

import argparse
import hashlib
import json
import math
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from nabu_protein.higher_order import additive_score
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

VERSION = "NABU_V9_1_DEGREE_NORMALIZED_PAIR_TRANSFER_DEV_V1"


def rank_hash(candidate_ids, prediction) -> str:
    frame = pd.DataFrame(
        {
            "candidate_id": list(candidate_ids),
            "prediction": np.asarray(prediction, dtype=float),
        }
    )
    order = (
        frame.sort_values(
            ["prediction", "candidate_id"],
            ascending=[False, True],
            kind="stable",
        )["candidate_id"]
        .astype(str)
        .tolist()
    )
    return hashlib.sha256("\n".join(order).encode("utf-8")).hexdigest()


def transferred_pair_effect(pair, factor):
    left, right = pair
    propensity = factor["propensity"]
    confidence = factor["confidence"]
    if left not in propensity or right not in propensity:
        return None
    raw = float(factor["bias"] + propensity[left] + propensity[right])
    weight = math.sqrt(float(confidence[left]) * float(confidence[right]))
    return float(raw * weight)


def score_v91(mutations, model, factor):
    base = model.model["base"]
    main_memory = base["main"]
    pair_memory = base["pair"]
    global_mean = float(base["global_mean"])

    if not all(mutation in main_memory for mutation in mutations):
        return global_mean, {
            "all_main_supported": False,
            "exact_pair_count": 0,
            "transferred_pair_count": 0,
            "unsupported_pair_count": (
                int(math.comb(len(mutations), 2))
                if len(mutations) >= 2
                else 0
            ),
            "transferred_active_node_count": 0,
            "transferred_mean_degree": 0.0,
            "raw_transferred_sum": 0.0,
            "normalized_transferred_sum": 0.0,
        }

    score = float(additive_score(mutations, global_mean, main_memory))
    exact_sum = 0.0
    transferred_values = []
    transferred_nodes = set()
    exact_count = 0
    unsupported_count = 0

    for pair in combinations(mutations, 2):
        exact = pair_memory.get(pair)
        if exact is not None:
            exact_sum += float(exact["mean"] * exact["confidence"])
            exact_count += 1
            continue

        transferred = transferred_pair_effect(pair, factor)
        if transferred is None:
            unsupported_count += 1
        else:
            transferred_values.append(float(transferred))
            transferred_nodes.update(pair)

    transferred_count = len(transferred_values)
    active_node_count = len(transferred_nodes)
    raw_sum = float(np.sum(transferred_values)) if transferred_values else 0.0

    if transferred_count > 0:
        mean_degree = float(
            (2.0 * transferred_count) / active_node_count
        )
        normalized_sum = float(raw_sum / mean_degree)
    else:
        mean_degree = 0.0
        normalized_sum = 0.0

    score = float(score + exact_sum + normalized_sum)

    return score, {
        "all_main_supported": True,
        "exact_pair_count": int(exact_count),
        "transferred_pair_count": int(transferred_count),
        "unsupported_pair_count": int(unsupported_count),
        "transferred_active_node_count": int(active_node_count),
        "transferred_mean_degree": float(mean_degree),
        "raw_transferred_sum": float(raw_sum),
        "normalized_transferred_sum": float(normalized_sum),
    }


def finite_spearman(target, prediction):
    value = float(spearmanr(target, prediction).statistic)
    return value if np.isfinite(value) else None


def summarize_scale(values):
    values = np.asarray(values, dtype=float)
    return {
        "mean": float(np.mean(values)),
        "std": float(np.std(values)),
        "min": float(np.min(values)),
        "max": float(np.max(values)),
    }


def summarize_diag(diag_frame: pd.DataFrame) -> dict:
    transferred_mask = diag_frame["transferred_pair_count"] > 0
    if bool(transferred_mask.any()):
        degree = diag_frame.loc[
            transferred_mask,
            "transferred_mean_degree",
        ].to_numpy(dtype=float)
        raw = np.abs(
            diag_frame.loc[
                transferred_mask,
                "raw_transferred_sum",
            ].to_numpy(dtype=float)
        )
        normalized = np.abs(
            diag_frame.loc[
                transferred_mask,
                "normalized_transferred_sum",
            ].to_numpy(dtype=float)
        )
        transfer_stats = {
            "mean_degree_mean": float(np.mean(degree)),
            "mean_degree_median": float(np.median(degree)),
            "mean_degree_max": float(np.max(degree)),
            "abs_raw_transfer_mean": float(np.mean(raw)),
            "abs_normalized_transfer_mean": float(np.mean(normalized)),
            "abs_raw_transfer_max": float(np.max(raw)),
            "abs_normalized_transfer_max": float(np.max(normalized)),
        }
    else:
        transfer_stats = {
            "mean_degree_mean": 0.0,
            "mean_degree_median": 0.0,
            "mean_degree_max": 0.0,
            "abs_raw_transfer_mean": 0.0,
            "abs_normalized_transfer_mean": 0.0,
            "abs_raw_transfer_max": 0.0,
            "abs_normalized_transfer_max": 0.0,
        }

    return {
        "all_main_supported_rows": int(
            diag_frame["all_main_supported"].sum()
        ),
        "rows_with_exact_pairs": int(
            (diag_frame["exact_pair_count"] > 0).sum()
        ),
        "rows_with_transferred_pairs": int(
            (diag_frame["transferred_pair_count"] > 0).sum()
        ),
        "exact_pair_contributions_total": int(
            diag_frame["exact_pair_count"].sum()
        ),
        "transferred_pair_contributions_total": int(
            diag_frame["transferred_pair_count"].sum()
        ),
        "unsupported_pair_contributions_total": int(
            diag_frame["unsupported_pair_count"].sum()
        ),
        **transfer_stats,
    }


def run(source_gz: Path, output_dir: Path) -> dict:
    source_sha, fit_frame, validation_frame, test_frame = load_ired(source_gz)
    reference = derive_reference(fit_frame["sequence"].tolist())

    fit_mutations = [
        mutation_set(sequence, reference)
        for sequence in fit_frame["sequence"]
    ]
    validation_mutations = [
        mutation_set(sequence, reference)
        for sequence in validation_frame["sequence"]
    ]
    test_mutations = [
        mutation_set(sequence, reference)
        for sequence in test_frame["sequence"]
    ]

    fit_target = fit_frame["target"].to_numpy(dtype=float)
    validation_target = validation_frame["target"].to_numpy(dtype=float)
    test_target = test_frame["target"].to_numpy(dtype=float)

    model = NabuV83Model().fit(
        mutation_sets=fit_mutations,
        labels=fit_target,
        candidate_ids=fit_frame["sequence"].tolist(),
    )

    factor = fit_pair_factorization(
        model,
        fit_mutations,
        fit_target,
        permute_targets=False,
    )
    permuted_factor = fit_pair_factorization(
        model,
        fit_mutations,
        fit_target,
        permute_targets=True,
    )

    test_ids = test_frame["sequence"].tolist()
    phase2_prediction, strict_scoreable = phase2_abstention_predictions(
        model,
        test_mutations,
        test_ids,
    )
    b2_prediction, all_main = b2_predictions(model, test_mutations)

    v90_rows = [
        score_v9(mutations, model, factor)
        for mutations in test_mutations
    ]
    v90_prediction = np.asarray(
        [value for value, _ in v90_rows],
        dtype=float,
    )

    v91_rows = [
        score_v91(mutations, model, factor)
        for mutations in test_mutations
    ]
    v91_prediction = np.asarray(
        [value for value, _ in v91_rows],
        dtype=float,
    )
    v91_diag = pd.DataFrame([diag for _, diag in v91_rows])

    permuted_rows = [
        score_v91(mutations, model, permuted_factor)
        for mutations in test_mutations
    ]
    permuted_prediction = np.asarray(
        [value for value, _ in permuted_rows],
        dtype=float,
    )

    replay_factor = fit_pair_factorization(
        model,
        fit_mutations,
        fit_target,
        permute_targets=False,
    )
    replay_prediction = np.asarray(
        [
            score_v91(mutations, model, replay_factor)[0]
            for mutations in test_mutations
        ],
        dtype=float,
    )

    arms = {
        "PHASE2_V8_3_ABSTENTION": phase2_prediction,
        "B2_MAIN_ONLY": b2_prediction,
        "V9_0_RAW_PAIR_TRANSFER": v90_prediction,
        "V9_1_DEGREE_NORMALIZED_TRANSFER": v91_prediction,
        "V9_1_DEGREE_NORMALIZED_PERMUTED_CONTROL": permuted_prediction,
    }

    if not all(np.isfinite(values).all() for values in arms.values()):
        raise RuntimeError("A V9.1 arm produced non-finite predictions.")

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

    non_strict_mask = all_main & ~strict_scoreable
    non_strict_metrics = {
        "count": int(non_strict_mask.sum()),
        "B2_MAIN_ONLY": evaluate(
            test_target[non_strict_mask],
            b2_prediction[non_strict_mask],
        ),
        "V9_0_RAW_PAIR_TRANSFER": evaluate(
            test_target[non_strict_mask],
            v90_prediction[non_strict_mask],
        ),
        "V9_1_DEGREE_NORMALIZED_TRANSFER": evaluate(
            test_target[non_strict_mask],
            v91_prediction[non_strict_mask],
        ),
        "V9_1_DEGREE_NORMALIZED_PERMUTED_CONTROL": evaluate(
            test_target[non_strict_mask],
            permuted_prediction[non_strict_mask],
        ),
    }

    validation_b2, _ = b2_predictions(model, validation_mutations)
    validation_v90 = np.asarray(
        [
            score_v9(mutations, model, factor)[0]
            for mutations in validation_mutations
        ],
        dtype=float,
    )
    validation_v91 = np.asarray(
        [
            score_v91(mutations, model, factor)[0]
            for mutations in validation_mutations
        ],
        dtype=float,
    )

    validation_metrics = {
        "B2_MAIN_ONLY": evaluate(validation_target, validation_b2),
        "V9_0_RAW_PAIR_TRANSFER": evaluate(
            validation_target,
            validation_v90,
        ),
        "V9_1_DEGREE_NORMALIZED_TRANSFER": evaluate(
            validation_target,
            validation_v91,
        ),
        "role": "POST_PHASE2_DEVELOPMENT_DIAGNOSTIC_ONLY",
    }

    b2 = metrics["B2_MAIN_ONLY"]
    v90 = metrics["V9_0_RAW_PAIR_TRANSFER"]
    v91 = metrics["V9_1_DEGREE_NORMALIZED_TRANSFER"]
    perm = metrics["V9_1_DEGREE_NORMALIZED_PERMUTED_CONTROL"]

    non_strict_b2 = non_strict_metrics["B2_MAIN_ONLY"]
    non_strict_v91 = non_strict_metrics[
        "V9_1_DEGREE_NORMALIZED_TRANSFER"
    ]

    checks = {
        "spearman_beats_b2": bool(
            v91["spearman"] is not None
            and b2["spearman"] is not None
            and v91["spearman"] > b2["spearman"]
        ),
        "spearman_beats_v9_0": bool(
            v91["spearman"] is not None
            and v90["spearman"] is not None
            and v91["spearman"] > v90["spearman"]
        ),
        "spearman_beats_permuted_control": bool(
            v91["spearman"] is not None
            and perm["spearman"] is not None
            and v91["spearman"] > perm["spearman"]
        ),
        "elite_metric_improves_over_b2": bool(
            v91["top1_percent_hits"] > b2["top1_percent_hits"]
            or v91["top1_percent_enrichment"]
            > b2["top1_percent_enrichment"]
        ),
        "normalized_regret_not_worse_than_b2": bool(
            v91["normalized_regret_top1pct"]
            <= b2["normalized_regret_top1pct"]
        ),
        "non_strict_region_spearman_beats_b2": bool(
            non_strict_v91["spearman"] is not None
            and non_strict_b2["spearman"] is not None
            and non_strict_v91["spearman"]
            > non_strict_b2["spearman"]
        ),
        "finite_predictions": bool(np.isfinite(v91_prediction).all()),
        "deterministic_rank_replay": bool(
            rank_hash(test_ids, v91_prediction)
            == rank_hash(test_ids, replay_prediction)
        ),
    }

    predictions = pd.DataFrame(
        {
            "candidate_id": test_ids,
            "mutation_count": [len(x) for x in test_mutations],
            "target": test_target,
            "strict_v83_scoreable": strict_scoreable,
            "all_main_supported": all_main,
            **arms,
            "v91_exact_pair_count": v91_diag["exact_pair_count"],
            "v91_transferred_pair_count": v91_diag[
                "transferred_pair_count"
            ],
            "v91_unsupported_pair_count": v91_diag[
                "unsupported_pair_count"
            ],
            "v91_transferred_active_node_count": v91_diag[
                "transferred_active_node_count"
            ],
            "v91_transferred_mean_degree": v91_diag[
                "transferred_mean_degree"
            ],
            "v91_raw_transferred_sum": v91_diag[
                "raw_transferred_sum"
            ],
            "v91_normalized_transferred_sum": v91_diag[
                "normalized_transferred_sum"
            ],
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
            "validation": int(len(validation_frame)),
            "test": int(len(test_frame)),
        },
        "reference_sha256": hashlib.sha256(
            reference.encode("utf-8")
        ).hexdigest(),
        "reference_length": int(len(reference)),
        "single_changed_mechanism": (
            "Transferred unseen-pair edge effects are divided by "
            "the mean degree of the transferred-edge subgraph."
        ),
        "factorization_reused_from_v9_0": True,
        "metrics": metrics,
        "non_strict_transfer_region_metrics": non_strict_metrics,
        "validation_diagnostic": validation_metrics,
        "by_mutation_count": by_mutation_count,
        "prediction_scale": {
            name: summarize_scale(values)
            for name, values in arms.items()
        },
        "transfer_diagnostics": summarize_diag(v91_diag),
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
            "V9_1_DEGREE_NORMALIZED_TRANSFER_PROMISING"
            if all(checks.values())
            else "V9_1_DEGREE_NORMALIZED_TRANSFER_NOT_YET_SUPPORTED"
        ),
        "interpretation_boundary": (
            "IRED is a post-blind development dataset. "
            "This result cannot repair Phase 2 and does not open Phase 3."
        ),
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "V9_1_RESULTS.json").write_text(
        json.dumps(result, indent=2, default=str),
        encoding="utf-8",
    )
    predictions.to_csv(
        output_dir / "V9_1_PREDICTIONS.csv",
        index=False,
    )

    print(json.dumps(result, indent=2, default=str))
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source_gz")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    run(Path(args.source_gz), Path(args.out))


if __name__ == "__main__":
    main()
