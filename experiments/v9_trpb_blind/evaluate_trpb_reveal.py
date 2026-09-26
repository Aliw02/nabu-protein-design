from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

DEV = Path(__file__).resolve().parent.parent / "v9_pair_transfer_dev"
sys.path.insert(0, str(DEV))

from run_v9_pair_transfer_ired import evaluate


PUBLISHED_RIDGE_SPEARMAN = 0.437
PUBLISHED_RIDGE_NDCG = 0.994
PUBLISHED_CARP_SUPERVISED_SPEARMAN = 0.509
PUBLISHED_CARP_SUPERVISED_NDCG = 0.996


def evaluate_by_count(target, prediction, counts):
    target = np.asarray(target, dtype=float)
    prediction = np.asarray(prediction, dtype=float)
    counts = np.asarray(counts, dtype=int)
    result = {}
    for count in sorted(set(counts.tolist())):
        mask = counts == count
        if int(mask.sum()) >= 2:
            result[str(count)] = evaluate(
                target[mask],
                prediction[mask],
            )
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("frozen_predictions_csv_gz")
    parser.add_argument("reveal_csv_gz")
    parser.add_argument("--freeze-manifest", required=True)
    parser.add_argument("--source-manifest", required=True)
    parser.add_argument("--train-diagnostics", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    pred = pd.read_csv(
        args.frozen_predictions_csv_gz,
        compression="gzip",
    )
    reveal = pd.read_csv(
        args.reveal_csv_gz,
        compression="gzip",
    )
    freeze = json.loads(
        Path(args.freeze_manifest).read_text(encoding="utf-8")
    )
    source = json.loads(
        Path(args.source_manifest).read_text(encoding="utf-8")
    )
    train = json.loads(
        Path(args.train_diagnostics).read_text(encoding="utf-8")
    )

    merged = pred.merge(
        reveal,
        on="row_id",
        how="inner",
        validate="one_to_one",
    )
    if len(merged) != source["observed_counts"]["test"]:
        raise RuntimeError("Reveal merge count mismatch.")

    target = merged["target"].to_numpy(dtype=float)
    counts = merged["mutation_count"].to_numpy(dtype=int)

    arm_names = [
        "B2_MAIN_ONLY",
        "V9_CONTEXTUAL_ESM_RESIDUAL",
        "V9_CONSENSUS_B2_CONTEXT",
    ]
    metrics = {
        arm: evaluate(
            target,
            merged[arm].to_numpy(dtype=float),
        )
        for arm in arm_names
    }
    by_count = {
        arm: evaluate_by_count(
            target,
            merged[arm].to_numpy(dtype=float),
            counts,
        )
        for arm in arm_names
    }

    clean_mask = merged["complete_clean_main"].astype(bool).to_numpy()
    if int(clean_mask.sum()) < 2:
        raise RuntimeError("Insufficient clean-main supported test rows.")

    clean_metrics = {
        arm: evaluate(
            target[clean_mask],
            merged.loc[clean_mask, arm].to_numpy(dtype=float),
        )
        for arm in arm_names
    }

    consensus = metrics["V9_CONSENSUS_B2_CONTEXT"]
    b2 = metrics["B2_MAIN_ONLY"]
    clean_consensus = clean_metrics["V9_CONSENSUS_B2_CONTEXT"]
    clean_b2 = clean_metrics["B2_MAIN_ONLY"]

    source_contract = bool(
        source["observed_counts"]
        == source["expected_counts"]
        and source["test_targets_excluded_from_prediction_artifact"]
        and freeze["source_sha256"] == source["source_sha256"]
        and freeze["test_count"] == source["observed_counts"]["test"]
        and freeze["freeze_completed_before_reveal"]
    )

    def positive_order(order: int) -> bool:
        item = by_count["V9_CONSENSUS_B2_CONTEXT"].get(str(order))
        return bool(
            item is not None
            and item["spearman"] is not None
            and item["spearman"] > 0.0
        )

    internal_checks = {
        "consensus_spearman_beats_b2": bool(
            consensus["spearman"] > b2["spearman"]
        ),
        "clean_region_spearman_beats_b2": bool(
            clean_consensus["spearman"] > clean_b2["spearman"]
        ),
        "top1_hits_at_least_b2": bool(
            consensus["top1_percent_hits"]
            >= b2["top1_percent_hits"]
        ),
        "regret_not_worse_than_b2": bool(
            consensus["normalized_regret_top1pct"]
            <= b2["normalized_regret_top1pct"]
        ),
        "positive_spearman_order_3": positive_order(3),
        "positive_spearman_order_4": positive_order(4),
        "frozen_hashes_present": bool(
            freeze.get("prediction_semantic_sha256")
            and freeze.get("b2_rank_sha256")
            and freeze.get("context_rank_sha256")
            and freeze.get("consensus_rank_sha256")
        ),
        "source_and_split_contracts_pass": source_contract,
        "validation_targets_unused_for_fit": bool(
            not train["validation_targets_used_for_fit"]
            and not source["validation_targets_used_for_fit"]
        ),
    }
    internal_pass = all(internal_checks.values())

    published_checks = {
        "internal_external_transfer_gate_pass": bool(internal_pass),
        "consensus_spearman_at_least_published_ridge": bool(
            consensus["spearman"] >= PUBLISHED_RIDGE_SPEARMAN
        ),
    }
    published_pass = all(published_checks.values())
    strong_published = bool(
        consensus["spearman"] >= PUBLISHED_CARP_SUPERVISED_SPEARMAN
    )

    if published_pass:
        decision = "TRPB_BLIND_PUBLISHED_BENCHMARK_PASS"
    elif internal_pass:
        decision = "TRPB_BLIND_INTERNAL_TRANSFER_PASS_BENCHMARK_FAIL"
    else:
        decision = "TRPB_BLIND_EXTERNAL_TRANSFER_FAIL"

    result = {
        "version": "NABU_V9_TRPB_BLIND_EVALUATION_V1",
        "scientific_claim": False,
        "phase3_opened": False,
        "dataset": "FLIP2_TrpB_two_to_many",
        "source_manifest": source,
        "freeze_manifest": freeze,
        "train_diagnostics": train,
        "coverage": {
            "test_count": int(len(merged)),
            "complete_clean_main": int(clean_mask.sum()),
            "b2_fallback_rows": int((~clean_mask).sum()),
        },
        "metrics": metrics,
        "complete_clean_main_region_metrics": clean_metrics,
        "by_mutation_count": by_count,
        "published_reference": {
            "Ridge_one_hot": {
                "spearman": PUBLISHED_RIDGE_SPEARMAN,
                "ndcg": PUBLISHED_RIDGE_NDCG,
            },
            "CARP_640M_supervised": {
                "spearman": PUBLISHED_CARP_SUPERVISED_SPEARMAN,
                "ndcg": PUBLISHED_CARP_SUPERVISED_NDCG,
            },
        },
        "internal_external_transfer_checks": internal_checks,
        "internal_external_transfer_pass": bool(internal_pass),
        "published_benchmark_checks": published_checks,
        "published_benchmark_pass": bool(published_pass),
        "strong_published_spearman_reference_reached": strong_published,
        "decision": decision,
        "boundary": (
            "This is the first TrpB reveal for this frozen architecture. "
            "No TrpB test target may be used for retuning after this result."
        ),
    }

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "TRPB_BLIND_RESULTS.json").write_text(
        json.dumps(result, indent=2),
        encoding="utf-8",
    )
    merged.to_csv(
        out / "TRPB_REVEALED_PREDICTIONS.csv.gz",
        index=False,
        compression="gzip",
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
