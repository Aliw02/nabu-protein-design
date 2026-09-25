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
PHASE2A2 = PHASE2A / "2A2_micro_batch_loop"
sys.path.insert(0, str(PHASE2A0))
sys.path.insert(0, str(PHASE2A1))
sys.path.insert(0, str(PHASE2A2))

from campaign import (  # noqa: E402
    CampaignSimulator,
    VirtualAssayOracle,
    identity_pool_from_frame,
    truth_frame_from_frame,
)
from baselines import HistoricalFiftyFiftyPolicy  # noqa: E402
from run_baseline_tournament import CampaignEvaluator  # noqa: E402
from run_micro_batch_ablation import fixed_batch_schedule  # noqa: E402
from adaptive_controller import EvidenceAdaptivePolicyV1  # noqa: E402


def json_default(value):
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    raise TypeError(f"Not JSON serializable: {type(value)!r}")


def build_policy(name: str, identity_pool: pd.DataFrame):
    if name == "historical_50_50":
        return HistoricalFiftyFiftyPolicy(identity_pool)
    if name == "evidence_adaptive_v1":
        return EvidenceAdaptivePolicyV1(identity_pool)
    raise ValueError(f"Unknown policy: {name}")


def run_condition(
    policy_name: str,
    batch_size: int,
    identity_pool: pd.DataFrame,
    truth: pd.DataFrame,
    output_dir: Path,
    seed_count: int,
    target_count: int,
) -> dict:
    condition_name = f"{policy_name}_micro_{batch_size}"
    condition_dir = output_dir / condition_name

    simulator = CampaignSimulator(
        identity_pool=identity_pool,
        oracle=VirtualAssayOracle(truth),
        output_dir=condition_dir / "campaign",
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
    policy = build_policy(policy_name, identity_pool)

    schedule = fixed_batch_schedule(
        seed_count=seed_count,
        target_count=target_count,
        batch_size=batch_size,
    )

    controller_rows = []

    for planned_batch in schedule:
        remaining = target_count - simulator.measurements_spent
        actual_batch = min(planned_batch, remaining)
        simulator.run_round(policy, batch_size=actual_batch)
        curve.append(evaluator.checkpoint(simulator))

        diagnostics = getattr(policy, "last_diagnostics", None)
        if diagnostics is not None:
            controller_rows.append(dict(diagnostics))

    curve_frame = pd.DataFrame(curve)
    curve_frame["condition"] = condition_name
    curve_frame.to_csv(
        condition_dir / "DISCOVERY_CURVE.csv",
        index=False,
    )

    if controller_rows:
        pd.DataFrame(controller_rows).to_csv(
            condition_dir / "CONTROLLER_TRAJECTORY.csv",
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
            "version": "NABU_PHASE2A3_CONTROLLER_RESULT_V1",
            "condition": condition_name,
            "policy": policy_name,
            "batch_size": int(batch_size),
            "seed_count": int(seed_count),
            "target_count": int(target_count),
            "refits_after_bootstrap": int(len(schedule)),
            "bootstrap_selection_sha256": bootstrap[
                "selection_sha256_before_label_reveal"
            ],
            "controller_trajectory_rows": int(len(controller_rows)),
            "scientific_scope": (
                "RhlA development-landscape controller ablation with "
                "cadence held fixed within each paired comparison."
            ),
        }
    )

    if controller_rows:
        pressures = np.asarray(
            [
                row["exploration_pressure"]
                for row in controller_rows
            ],
            dtype=float,
        )
        summary["controller_pressure"] = {
            "initial": float(pressures[0]),
            "final": float(pressures[-1]),
            "min": float(np.min(pressures)),
            "max": float(np.max(pressures)),
            "mean": float(np.mean(pressures)),
        }

    (condition_dir / "CONDITION_SUMMARY.json").write_text(
        json.dumps(summary, indent=2, default=json_default),
        encoding="utf-8",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run Phase-2A.3 evidence-adaptive controller ablation."
    )
    parser.add_argument(
        "--input",
        default="rhla_sample_efficiency_sealed_input/TRAINING_POOL.csv",
    )
    parser.add_argument(
        "--out",
        default=(
            "phase2/2A_closed_loop_acquisition/"
            "2A3_adaptive_controller/results/rhla_controller_v1"
        ),
    )
    parser.add_argument("--seed-fraction", type=float, default=0.05)
    parser.add_argument("--target-fraction", type=float, default=0.20)
    args = parser.parse_args()

    if not (0.0 < args.seed_fraction < args.target_fraction <= 1.0):
        raise SystemExit(
            "Require 0 < seed-fraction < target-fraction <= 1."
        )

    source = pd.read_csv(args.input)
    identity_pool = identity_pool_from_frame(source)
    truth = truth_frame_from_frame(source)

    pool_size = len(identity_pool)
    seed_count = int(math.ceil(pool_size * args.seed_fraction))
    target_count = int(math.ceil(pool_size * args.target_fraction))

    matrix = [
        ("historical_50_50", 8),
        ("evidence_adaptive_v1", 8),
        ("historical_50_50", 16),
        ("evidence_adaptive_v1", 16),
    ]

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    summaries = []
    for policy_name, batch_size in matrix:
        print(f"[2A3] START {policy_name} micro_{batch_size}")
        summary = run_condition(
            policy_name=policy_name,
            batch_size=batch_size,
            identity_pool=identity_pool,
            truth=truth,
            output_dir=out,
            seed_count=seed_count,
            target_count=target_count,
        )
        summaries.append(summary)
        print(
            f"[2A3] DONE {policy_name} micro_{batch_size} "
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
            "Controller conditions did not share the exact same bootstrap."
        )

    comparison_rows = []
    for summary in summaries:
        final = summary["final"]
        pressure = summary.get("controller_pressure", {})
        comparison_rows.append(
            {
                "condition": summary["condition"],
                "policy": summary["policy"],
                "batch_size": summary["batch_size"],
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
                "pressure_initial": pressure.get("initial"),
                "pressure_final": pressure.get("final"),
                "pressure_mean": pressure.get("mean"),
            }
        )

    comparison = pd.DataFrame(comparison_rows)
    comparison.to_csv(
        out / "CONTROLLER_COMPARISON.csv",
        index=False,
    )

    manifest = {
        "version": "NABU_PHASE2A3_CONTROLLER_ABLATION_V1",
        "scientific_claim": False,
        "stage": "2A3_EVIDENCE_ADAPTIVE_CONTROLLER_DEVELOPMENT",
        "controller_formula_frozen_before_run": True,
        "input": str(args.input),
        "pool_size": int(pool_size),
        "seed_fraction": float(args.seed_fraction),
        "target_fraction": float(args.target_fraction),
        "seed_count": int(seed_count),
        "target_count": int(target_count),
        "primary_cadence": 8,
        "sensitivity_cadence": 16,
        "shared_bootstrap_selection_sha256": next(iter(bootstrap_hashes)),
        "comparison": comparison_rows,
        "interpretation_rule": (
            "RhlA is a development landscape. Preserve failures. "
            "Do not claim general superiority or retune V1 after seeing "
            "these outcomes without versioning a new controller."
        ),
    }
    (out / "CONTROLLER_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, default=json_default),
        encoding="utf-8",
    )

    print(json.dumps(manifest, indent=2, default=json_default))


if __name__ == "__main__":
    main()
