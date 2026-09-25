from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
PHASE2A = HERE.parent
PHASE2A0 = PHASE2A / "2A0_campaign_simulator"
PHASE2A1 = PHASE2A / "2A1_baselines"
sys.path.insert(0, str(PHASE2A0))
sys.path.insert(0, str(PHASE2A1))

from campaign import (  # noqa: E402
    CampaignSimulator,
    VirtualAssayOracle,
    identity_pool_from_frame,
    truth_frame_from_frame,
)
from baselines import HistoricalFiftyFiftyPolicy  # noqa: E402
from run_baseline_tournament import CampaignEvaluator  # noqa: E402


def json_default(value):
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    raise TypeError(f"Not JSON serializable: {type(value)!r}")


def fixed_batch_schedule(
    seed_count: int,
    target_count: int,
    batch_size: int,
) -> list[int]:
    if seed_count >= target_count:
        raise ValueError("seed_count must be smaller than target_count.")
    if batch_size < 1:
        raise ValueError("batch_size must be positive.")

    schedule = []
    current = int(seed_count)
    while current < target_count:
        addition = min(int(batch_size), target_count - current)
        schedule.append(addition)
        current += addition
    return schedule


def historical_jump_schedule(
    seed_count: int,
    ten_percent_count: int,
    target_count: int,
) -> list[int]:
    if not (seed_count < ten_percent_count < target_count):
        raise ValueError(
            "Require seed_count < ten_percent_count < target_count."
        )
    return [
        int(ten_percent_count - seed_count),
        int(target_count - ten_percent_count),
    ]


def build_schedules(
    seed_count: int,
    ten_percent_count: int,
    target_count: int,
) -> dict[str, list[int]]:
    return {
        "historical_jumps": historical_jump_schedule(
            seed_count,
            ten_percent_count,
            target_count,
        ),
        "micro_32": fixed_batch_schedule(seed_count, target_count, 32),
        "micro_16": fixed_batch_schedule(seed_count, target_count, 16),
        "micro_8": fixed_batch_schedule(seed_count, target_count, 8),
        "micro_4": fixed_batch_schedule(seed_count, target_count, 4),
    }


def run_cadence(
    cadence_name: str,
    schedule: list[int],
    identity_pool: pd.DataFrame,
    truth: pd.DataFrame,
    output_dir: Path,
    seed_count: int,
    target_count: int,
) -> dict:
    cadence_dir = output_dir / cadence_name
    simulator = CampaignSimulator(
        identity_pool=identity_pool,
        oracle=VirtualAssayOracle(truth),
        output_dir=cadence_dir / "campaign",
    )
    bootstrap = simulator.bootstrap(
        seed_count=seed_count,
        min_per_fold=1,
    )
    bootstrap_ids = set(bootstrap["candidate_ids"])

    evaluator = CampaignEvaluator(
        identity_pool=identity_pool,
        truth=truth,
        bootstrap_count=simulator.measurements_spent,
    )
    curve = [evaluator.checkpoint(simulator)]

    policy = HistoricalFiftyFiftyPolicy(identity_pool)

    for batch_size in schedule:
        if simulator.measurements_spent >= target_count:
            break
        expected = min(
            int(batch_size),
            target_count - simulator.measurements_spent,
        )
        simulator.run_round(policy, batch_size=expected)
        curve.append(evaluator.checkpoint(simulator))

    if simulator.measurements_spent != target_count:
        raise RuntimeError(
            f"{cadence_name} stopped at {simulator.measurements_spent}, "
            f"expected {target_count}."
        )

    curve_frame = pd.DataFrame(curve)
    curve_frame["cadence"] = cadence_name
    curve_frame.to_csv(
        cadence_dir / "DISCOVERY_CURVE.csv",
        index=False,
    )

    final_ids = simulator.measured["candidate_id"].astype(str).tolist()
    summary = evaluator.summarize(
        curve=curve,
        bootstrap_ids=bootstrap_ids,
        final_measured_ids=final_ids,
    )
    summary.update(
        {
            "version": "NABU_PHASE2A2_CADENCE_RESULT_V1",
            "cadence": cadence_name,
            "policy": "historical_50_50",
            "seed_count": int(seed_count),
            "target_count": int(target_count),
            "batch_schedule": [int(value) for value in schedule],
            "refits_after_bootstrap": int(len(schedule)),
            "max_batch_size": int(max(schedule)),
            "min_batch_size": int(min(schedule)),
            "bootstrap_selection_sha256": bootstrap[
                "selection_sha256_before_label_reveal"
            ],
            "scientific_scope": (
                "RhlA development-landscape refit-cadence ablation under "
                "a fixed historical 50/50 acquisition policy."
            ),
        }
    )

    (cadence_dir / "CADENCE_SUMMARY.json").write_text(
        json.dumps(summary, indent=2, default=json_default),
        encoding="utf-8",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run Phase-2A.2 refit-cadence ablation with the acquisition "
            "policy held fixed."
        )
    )
    parser.add_argument(
        "--input",
        default="rhla_sample_efficiency_sealed_input/TRAINING_POOL.csv",
    )
    parser.add_argument(
        "--out",
        default=(
            "phase2/2A_closed_loop_acquisition/"
            "2A2_micro_batch_loop/results/rhla_micro_batch_ablation"
        ),
    )
    parser.add_argument("--seed-fraction", type=float, default=0.05)
    parser.add_argument("--middle-fraction", type=float, default=0.10)
    parser.add_argument("--target-fraction", type=float, default=0.20)
    args = parser.parse_args()

    if not (
        0.0
        < args.seed_fraction
        < args.middle_fraction
        < args.target_fraction
        <= 1.0
    ):
        raise SystemExit(
            "Require 0 < seed-fraction < middle-fraction "
            "< target-fraction <= 1."
        )

    source = pd.read_csv(args.input)
    identity_pool = identity_pool_from_frame(source)
    truth = truth_frame_from_frame(source)

    pool_size = len(identity_pool)
    seed_count = int(math.ceil(pool_size * args.seed_fraction))
    middle_count = int(math.ceil(pool_size * args.middle_fraction))
    target_count = int(math.ceil(pool_size * args.target_fraction))

    schedules = build_schedules(
        seed_count=seed_count,
        ten_percent_count=middle_count,
        target_count=target_count,
    )

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    summaries = []
    for cadence_name, schedule in schedules.items():
        print(
            f"[2A2] START {cadence_name} "
            f"refits={len(schedule)} max_batch={max(schedule)}"
        )
        summary = run_cadence(
            cadence_name=cadence_name,
            schedule=schedule,
            identity_pool=identity_pool,
            truth=truth,
            output_dir=out,
            seed_count=seed_count,
            target_count=target_count,
        )
        summaries.append(summary)
        print(
            f"[2A2] DONE {cadence_name} "
            f"best={summary['final']['best_true_fitness']:.6f} "
            f"top1={summary['final']['cumulative_top1_hits']} "
            f"auc={summary['post_bootstrap_normalized_discovery_auc']:.6f}"
        )

    bootstrap_hashes = {
        summary["bootstrap_selection_sha256"]
        for summary in summaries
    }
    if len(bootstrap_hashes) != 1:
        raise RuntimeError(
            "Cadence conditions did not share the exact same bootstrap."
        )

    comparison_rows = []
    for summary in summaries:
        final = summary["final"]
        comparison_rows.append(
            {
                "cadence": summary["cadence"],
                "policy": summary["policy"],
                "refits_after_bootstrap": summary[
                    "refits_after_bootstrap"
                ],
                "max_batch_size": summary["max_batch_size"],
                "final_best_true_fitness": final["best_true_fitness"],
                "final_regret": final["best_so_far_regret"],
                "final_top1_hits": final["cumulative_top1_hits"],
                "acquired_top1_hits_after_bootstrap": summary[
                    "acquired_top1_hits_after_bootstrap"
                ],
                "first_top1_measurement": summary[
                    "first_top1_measurement"
                ],
                "post_bootstrap_normalized_discovery_auc": summary[
                    "post_bootstrap_normalized_discovery_auc"
                ],
                "final_scoreable_coverage": final[
                    "scoreable_coverage_unmeasured"
                ],
                "final_router_mode": final["router_mode"],
                "final_visible_b3_oof_spearman": final[
                    "visible_b3_oof_spearman"
                ],
            }
        )

    comparison = pd.DataFrame(comparison_rows)
    comparison.to_csv(out / "MICRO_BATCH_COMPARISON.csv", index=False)

    manifest = {
        "version": "NABU_PHASE2A2_MICRO_BATCH_ABLATION_V1",
        "scientific_claim": False,
        "stage": "2A2_BATCH_STALENESS_CHARACTERIZATION",
        "policy_held_fixed": "historical_50_50",
        "input": str(args.input),
        "pool_size": int(pool_size),
        "seed_fraction": float(args.seed_fraction),
        "middle_fraction": float(args.middle_fraction),
        "target_fraction": float(args.target_fraction),
        "seed_count": int(seed_count),
        "middle_count": int(middle_count),
        "target_count": int(target_count),
        "shared_bootstrap_selection_sha256": next(iter(bootstrap_hashes)),
        "schedules": schedules,
        "comparison": comparison_rows,
        "interpretation_rule": (
            "Differences among rows isolate refit cadence under the same "
            "historical 50/50 acquisition rule. Do not infer a universal "
            "optimal batch size from RhlA alone."
        ),
    }
    (out / "MICRO_BATCH_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, default=json_default),
        encoding="utf-8",
    )

    print(json.dumps(manifest, indent=2, default=json_default))


if __name__ == "__main__":
    main()
