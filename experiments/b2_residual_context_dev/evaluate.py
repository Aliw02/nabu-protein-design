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
    parser.add_argument("--train-diagnostics", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    pred = pd.read_csv(args.frozen_predictions, compression="gzip")
    reveal = pd.read_csv(args.reveal, compression="gzip")
    freeze = json.loads(
        Path(args.freeze_manifest).read_text(encoding="utf-8")
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
    if len(frame) != freeze["test_count"]:
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

    def positive_order(order):
        item = by_count["B2_PLUS_CONTEXT_RESIDUAL"].get(str(order))
        return bool(
            item is not None
            and item["spearman"] is not None
            and item["spearman"] > 0
        )

    checks = {
        "full_spearman_beats_b2": bool(
            corrected["spearman"] > b2["spearman"]
        ),
        "top1_hits_at_least_b2": bool(
            corrected["top1_percent_hits"] >= b2["top1_percent_hits"]
        ),
        "regret_not_worse_than_b2": bool(
            corrected["normalized_regret_top1pct"]
            <= b2["normalized_regret_top1pct"]
        ),
        "positive_spearman_order_3": positive_order(3),
        "positive_spearman_order_4": positive_order(4),
        "beats_shuffled_context": bool(
            corrected["spearman"] > shuffled["spearman"]
        ),
        "freeze_before_reveal": bool(
            freeze["freeze_completed_before_reveal"]
        ),
    }
    passed = all(checks.values())

    result = {
        "version": "NABU_B2_RESIDUAL_CONTEXT_TRPB_DEV_V1",
        "status": "DEVELOPMENT_AFTER_TRPB_REVEAL",
        "scientific_claim": False,
        "phase3_opened": False,
        "train_diagnostics": train,
        "freeze_manifest": freeze,
        "metrics": metrics,
        "by_mutation_count": by_count,
        "decision_checks": checks,
        "development_pass": bool(passed),
        "decision": (
            "B2_RESIDUAL_CONTEXT_ELIGIBLE_FOR_UNTOUCHED_FREEZE"
            if passed
            else "B2_RESIDUAL_CONTEXT_REJECT"
        ),
    }

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "RESULTS.json").write_text(
        json.dumps(result, indent=2),
        encoding="utf-8",
    )
    frame.to_csv(
        out / "REVEALED_PREDICTIONS.csv.gz",
        index=False,
        compression="gzip",
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
