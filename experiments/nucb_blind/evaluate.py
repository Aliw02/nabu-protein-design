from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

DEV = Path(__file__).resolve().parent.parent / "v9_pair_transfer_dev"
sys.path.insert(0, str(DEV))

from run_v9_pair_transfer_ired import evaluate


PUBLISHED_RIDGE_SPEARMAN = 0.606
PUBLISHED_RIDGE_NDCG = 0.943
PUBLISHED_RIDGE_LIKELIHOOD_SPEARMAN = 0.623
PUBLISHED_RIDGE_LIKELIHOOD_NDCG = 0.946
PUBLISHED_CARP_SUPERVISED_SPEARMAN = 0.717
PUBLISHED_CARP_SUPERVISED_NDCG = 0.970
PUBLISHED_ESMC_SUPERVISED_SPEARMAN = 0.723
PUBLISHED_ESMC_SUPERVISED_NDCG = 0.966


def evaluate_by_count(target, prediction, counts):
    result = {}
    counts = np.asarray(counts, dtype=int)
    for count in sorted(set(counts.tolist())):
        mask = counts == count
        if int(mask.sum()) >= 2:
            result[str(count)] = evaluate(
                np.asarray(target)[mask],
                np.asarray(prediction)[mask],
            )
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("frozen_predictions")
    parser.add_argument("reveal")
    parser.add_argument("--freeze-manifest", required=True)
    parser.add_argument("--replay-freeze-manifest", required=True)
    parser.add_argument("--source-manifest", required=True)
    parser.add_argument("--train-diagnostics", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    pred = pd.read_csv(args.frozen_predictions, compression="gzip")
    reveal = pd.read_csv(args.reveal, compression="gzip")
    freeze = json.loads(
        Path(args.freeze_manifest).read_text(encoding="utf-8")
    )
    replay = json.loads(
        Path(args.replay_freeze_manifest).read_text(encoding="utf-8")
    )
    source = json.loads(
        Path(args.source_manifest).read_text(encoding="utf-8")
    )
    train = json.loads(
        Path(args.train_diagnostics).read_text(encoding="utf-8")
    )

    frame = pred.merge(
        reveal,
        on="row_id",
        how="inner",
        validate="one_to_one",
    )
    if len(frame) != int(source["observed_counts"]["test"]):
        raise RuntimeError("Reveal merge count mismatch.")

    target = frame["target"].to_numpy(dtype=float)
    counts = frame["mutation_count"].to_numpy(dtype=int)

    arms = [
        "B2_MAIN_ONLY",
        "B2_PLUS_CONTEXT_RESIDUAL",
        "B2_PLUS_SHUFFLED_CONTEXT",
    ]
    metrics = {
        arm: evaluate(target, frame[arm].to_numpy(dtype=float))
        for arm in arms
    }
    by_count = {
        arm: evaluate_by_count(
            target,
            frame[arm].to_numpy(dtype=float),
            counts,
        )
        for arm in arms
    }

    b2 = metrics["B2_MAIN_ONLY"]
    corrected = metrics["B2_PLUS_CONTEXT_RESIDUAL"]
    shuffled = metrics["B2_PLUS_SHUFFLED_CONTEXT"]

    relevant_order_checks = {}
    for order, item in by_count["B2_PLUS_CONTEXT_RESIDUAL"].items():
        if int(item["count"]) >= 100:
            relevant_order_checks[order] = bool(
                item["spearman"] is not None
                and item["spearman"] > 0.0
            )

    determinism = bool(
        freeze["prediction_semantic_sha256"]
        == replay["prediction_semantic_sha256"]
    )

    source_contract = bool(
        source["source_md5"] == "cdedaadd57f2b148950b357e9f24c238"
        and source["test_targets_excluded_from_prediction_artifact"]
        and source["scope_contract_pass"]
        and all(
            int(order) <= 2
            for order in source["fit_mutation_order_distribution"]
        )
        and all(
            int(order) >= 3
            for order in source["test_mutation_order_distribution"]
        )
    )

    finite_predictions = bool(
        np.isfinite(
            frame[
                [
                    "B2_MAIN_ONLY",
                    "B2_PLUS_CONTEXT_RESIDUAL",
                    "B2_PLUS_SHUFFLED_CONTEXT",
                ]
            ].to_numpy(dtype=float)
        ).all()
    )

    checks = {
        "corrected_spearman_beats_b2": bool(
            corrected["spearman"] > b2["spearman"]
        ),
        "corrected_ndcg_at_least_b2": bool(
            corrected["ndcg"] >= b2["ndcg"]
        ),
        "corrected_spearman_beats_shuffled": bool(
            corrected["spearman"] > shuffled["spearman"]
        ),
        "published_ridge_spearman_reached": bool(
            corrected["spearman"] >= PUBLISHED_RIDGE_SPEARMAN
        ),
        "freeze_before_reveal": bool(
            freeze["freeze_completed_before_reveal"]
        ),
        "source_and_split_semantics_pass": source_contract,
        "validation_targets_unused_for_fit": bool(
            not source["validation_targets_used_for_fit"]
            and not train["validation_targets_used_for_fit"]
        ),
        "finite_predictions": finite_predictions,
        "deterministic_replay": determinism,
        "positive_spearman_all_orders_n_ge_100": bool(
            relevant_order_checks
            and all(relevant_order_checks.values())
        ),
    }
    passed = all(checks.values())

    result = {
        "version": "NABU_B2_RESIDUAL_CONTEXT_NUCB_BLIND_V1",
        "dataset": "FLIP2_NucB_two_to_many",
        "scientific_claim": False,
        "phase3_opened": False,
        "source_manifest": source,
        "train_diagnostics": train,
        "freeze_manifest": freeze,
        "replay_freeze_manifest": replay,
        "metrics": metrics,
        "by_mutation_count": by_count,
        "positive_order_checks_n_ge_100": relevant_order_checks,
        "published_reference": {
            "Ridge_one_hot": {
                "spearman": PUBLISHED_RIDGE_SPEARMAN,
                "ndcg": PUBLISHED_RIDGE_NDCG,
            },
            "Ridge_one_hot_plus_likelihoods": {
                "spearman": PUBLISHED_RIDGE_LIKELIHOOD_SPEARMAN,
                "ndcg": PUBLISHED_RIDGE_LIKELIHOOD_NDCG,
            },
            "CARP_640M_supervised": {
                "spearman": PUBLISHED_CARP_SUPERVISED_SPEARMAN,
                "ndcg": PUBLISHED_CARP_SUPERVISED_NDCG,
            },
            "ESMC_300M_supervised": {
                "spearman": PUBLISHED_ESMC_SUPERVISED_SPEARMAN,
                "ndcg": PUBLISHED_ESMC_SUPERVISED_NDCG,
            },
        },
        "decision_checks": checks,
        "blind_pass": bool(passed),
        "strong_reference_reached": bool(
            corrected["spearman"] >= PUBLISHED_ESMC_SUPERVISED_SPEARMAN
        ),
        "decision": (
            "NUCB_UNTOUCHED_EXTERNAL_PASS"
            if passed
            else "NUCB_UNTOUCHED_EXTERNAL_FAIL"
        ),
        "boundary": (
            "This is the first NucB reveal in the project. "
            "No NucB target may be used for post-reveal retuning."
        ),
    }

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "NUCB_BLIND_RESULTS.json").write_text(
        json.dumps(result, indent=2),
        encoding="utf-8",
    )
    frame.to_csv(
        out / "NUCB_REVEALED_PREDICTIONS.csv.gz",
        index=False,
        compression="gzip",
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
