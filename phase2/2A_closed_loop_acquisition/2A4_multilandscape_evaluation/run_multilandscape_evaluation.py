from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
PHASE2A = HERE.parent
PHASE2A0 = PHASE2A / "2A0_campaign_simulator"
PHASE2A1 = PHASE2A / "2A1_baselines"
PHASE2A2 = PHASE2A / "2A2_micro_batch_loop"
PHASE2A3 = PHASE2A / "2A3_adaptive_controller"

for path in (PHASE2A0, PHASE2A1, PHASE2A2, PHASE2A3):
    sys.path.insert(0, str(path))

from campaign import (  # noqa: E402
    CampaignSimulator,
    VirtualAssayOracle,
    identity_pool_from_frame,
    truth_frame_from_frame,
)
from baselines import (  # noqa: E402
    DeterministicRandomPolicy,
    HistoricalFiftyFiftyPolicy,
)
from run_baseline_tournament import CampaignEvaluator  # noqa: E402
from run_micro_batch_ablation import fixed_batch_schedule  # noqa: E402
from adaptive_controller import RelativeEvidencePolicyV2  # noqa: E402
from nabu_protein.higher_order import fnv1a32  # noqa: E402


POOL_NAMESPACE = "NABU_PHASE2A4_IDENTITY_POOL_V1"


def json_default(value):
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    raise TypeError(f"Not JSON serializable: {type(value)!r}")


def load_standard_landscape(path: str | Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    required = {"mutant", "DMS_score"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(
            f"Landscape {path} is missing columns: {sorted(missing)}"
        )

    data = frame[["mutant", "DMS_score"]].dropna().copy()
    data["mutant"] = data["mutant"].astype(str)
    data["DMS_score"] = pd.to_numeric(
        data["DMS_score"],
        errors="coerce",
    )
    data = data.dropna(subset=["DMS_score"]).copy()
    data["candidate_id"] = data["mutant"].astype(str)

    if data["candidate_id"].duplicated().any():
        data = (
            data.groupby("candidate_id", as_index=False, sort=False)
            .agg({"mutant": "first", "DMS_score": "mean"})
        )

    return data.reset_index(drop=True)


def deterministic_identity_subpool(
    frame: pd.DataFrame,
    pool_size: int,
) -> pd.DataFrame:
    if pool_size < 1:
        raise ValueError("pool_size must be positive.")
    if len(frame) < pool_size:
        raise ValueError(
            f"Source landscape has {len(frame)} candidates, "
            f"smaller than requested pool {pool_size}."
        )

    work = frame.copy()
    work["_pool_hash"] = work["candidate_id"].astype(str).map(
        lambda candidate_id: fnv1a32(
            f"{POOL_NAMESPACE}|{candidate_id}"
        )
    )
    work = work.sort_values(
        ["_pool_hash", "candidate_id"],
        ascending=[True, True],
    ).head(pool_size)
    return work.drop(columns=["_pool_hash"]).reset_index(drop=True)


def build_policy(name: str, identity_pool: pd.DataFrame):
    if name == "deterministic_random":
        return DeterministicRandomPolicy(seed=161)
    if name == "historical_50_50":
        return HistoricalFiftyFiftyPolicy(identity_pool)
    if name == "relative_evidence_v2":
        return RelativeEvidencePolicyV2(identity_pool)
    raise ValueError(f"Unknown policy: {name}")


def run_condition(
    landscape_name: str,
    policy_name: str,
    batch_size: int,
    identity_pool: pd.DataFrame,
    truth: pd.DataFrame,
    output_dir: Path,
    seed_count: int,
    target_count: int,
) -> dict:
    condition_name = f"{policy_name}_micro_{batch_size}"
    condition_dir = output_dir / landscape_name / condition_name

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
    curve_frame["landscape"] = landscape_name
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
            "version": "NABU_PHASE2A4_CONDITION_RESULT_V1",
            "landscape": landscape_name,
            "condition": condition_name,
            "policy": policy_name,
            "batch_size": int(batch_size),
            "seed_count": int(seed_count),
            "target_count": int(target_count),
            "bootstrap_selection_sha256": bootstrap[
                "selection_sha256_before_label_reveal"
            ],
            "rounds_completed": int(len(schedule)),
        }
    )

    if controller_rows:
        pressures = np.asarray(
            [row["exploration_pressure"] for row in controller_rows],
            dtype=float,
        )
        summary["controller_pressure"] = {
            "initial": float(pressures[0]),
            "final": float(pressures[-1]),
            "mean": float(np.mean(pressures)),
            "min": float(np.min(pressures)),
            "max": float(np.max(pressures)),
        }

    (condition_dir / "CONDITION_SUMMARY.json").write_text(
        json.dumps(summary, indent=2, default=json_default),
        encoding="utf-8",
    )
    return summary


def parse_landscape(value: str) -> tuple[str, str]:
    if "=" not in value:
        raise argparse.ArgumentTypeError(
            "Landscape must be supplied as NAME=PATH."
        )
    name, path = value.split("=", 1)
    name = name.strip()
    path = path.strip()
    if not name or not path:
        raise argparse.ArgumentTypeError(
            "Landscape must be supplied as NAME=PATH."
        )
    return name, path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run preregistered NABU Phase-2A.4 evaluation."
    )
    parser.add_argument(
        "--landscape",
        action="append",
        required=True,
        type=parse_landscape,
        help="Repeat as NAME=PATH.",
    )
    parser.add_argument(
        "--out",
        default=(
            "phase2/2A_closed_loop_acquisition/"
            "2A4_multilandscape_evaluation/results/"
            "multilandscape_v1"
        ),
    )
    parser.add_argument("--pool-size", type=int, default=2048)
    parser.add_argument("--seed-count", type=int, default=64)
    parser.add_argument("--target-count", type=int, default=256)
    args = parser.parse_args()

    if not (0 < args.seed_count < args.target_count <= args.pool_size):
        raise SystemExit(
            "Require 0 < seed-count < target-count <= pool-size."
        )

    matrix = [
        ("deterministic_random", 8),
        ("historical_50_50", 8),
        ("relative_evidence_v2", 8),
        ("historical_50_50", 16),
        ("relative_evidence_v2", 16),
    ]

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    summaries = []
    pool_manifests = []

    for landscape_name, landscape_path in args.landscape:
        source = load_standard_landscape(landscape_path)
        subpool = deterministic_identity_subpool(
            source,
            pool_size=args.pool_size,
        )

        identity_pool = identity_pool_from_frame(subpool)
        truth = truth_frame_from_frame(subpool)

        pool_manifest = {
            "landscape": landscape_name,
            "source_rows": int(len(source)),
            "campaign_pool_rows": int(len(subpool)),
            "pool_namespace": POOL_NAMESPACE,
            "candidate_ids_sha256_basis": (
                "sorted by FNV-1a identity hash then candidate_id"
            ),
        }
        pool_manifests.append(pool_manifest)

        for policy_name, batch_size in matrix:
            print(
                f"[2A4] START {landscape_name} "
                f"{policy_name} micro_{batch_size}"
            )
            summary = run_condition(
                landscape_name=landscape_name,
                policy_name=policy_name,
                batch_size=batch_size,
                identity_pool=identity_pool,
                truth=truth,
                output_dir=out,
                seed_count=args.seed_count,
                target_count=args.target_count,
            )
            summaries.append(summary)
            print(
                f"[2A4] DONE {landscape_name} "
                f"{policy_name} micro_{batch_size} "
                f"best={summary['final']['best_true_fitness']:.6f} "
                f"top1={summary['final']['cumulative_top1_hits']} "
                f"auc={summary['post_bootstrap_normalized_discovery_auc']:.6f}"
            )

    rows = []
    for summary in summaries:
        final = summary["final"]
        pressure = summary.get("controller_pressure", {})
        rows.append(
            {
                "landscape": summary["landscape"],
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

    comparison = pd.DataFrame(rows)
    comparison.to_csv(
        out / "MULTILANDSCAPE_COMPARISON.csv",
        index=False,
    )

    paired_rows = []
    for landscape_name, _ in args.landscape:
        for batch_size in (8, 16):
            subset = comparison[
                (comparison["landscape"] == landscape_name)
                & (comparison["batch_size"] == batch_size)
                & comparison["policy"].isin(
                    ["historical_50_50", "relative_evidence_v2"]
                )
            ].set_index("policy")

            if len(subset) != 2:
                raise RuntimeError(
                    f"Missing paired controller rows for "
                    f"{landscape_name}/micro_{batch_size}."
                )

            reference = subset.loc["historical_50_50"]
            adaptive = subset.loc["relative_evidence_v2"]
            paired_rows.append(
                {
                    "landscape": landscape_name,
                    "batch_size": int(batch_size),
                    "delta_auc_v2_minus_50_50": float(
                        adaptive[
                            "post_bootstrap_normalized_discovery_auc"
                        ]
                        - reference[
                            "post_bootstrap_normalized_discovery_auc"
                        ]
                    ),
                    "delta_top1_hits_v2_minus_50_50": int(
                        adaptive["final_top1_hits"]
                        - reference["final_top1_hits"]
                    ),
                    "delta_regret_v2_minus_50_50": float(
                        adaptive["final_regret"]
                        - reference["final_regret"]
                    ),
                    "delta_first_top1_measurement_v2_minus_50_50": (
                        None
                        if pd.isna(adaptive["first_top1_measurement"])
                        or pd.isna(reference["first_top1_measurement"])
                        else int(
                            adaptive["first_top1_measurement"]
                            - reference["first_top1_measurement"]
                        )
                    ),
                }
            )

    paired = pd.DataFrame(paired_rows)
    paired.to_csv(
        out / "PAIRED_CONTROLLER_DELTAS.csv",
        index=False,
    )

    manifest = {
        "version": "NABU_PHASE2A4_MULTILANDSCAPE_EVALUATION_V1",
        "scientific_claim": False,
        "stage": "2A4_MULTILANDSCAPE_GENERALIZATION",
        "controllers_frozen_before_run": True,
        "per_landscape_retuning": False,
        "landscapes": [name for name, _ in args.landscape],
        "pool_size": int(args.pool_size),
        "seed_count": int(args.seed_count),
        "target_count": int(args.target_count),
        "primary_cadence": 8,
        "sensitivity_cadence": 16,
        "pool_manifests": pool_manifests,
        "comparison": rows,
        "paired_controller_deltas": paired_rows,
        "interpretation_rule": (
            "Evaluate paired controller behavior across landscapes. "
            "These historically used landscapes are not the final "
            "untouched Phase-2D transfer set."
        ),
    }
    (out / "MULTILANDSCAPE_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, default=json_default),
        encoding="utf-8",
    )

    print(json.dumps(manifest, indent=2, default=json_default))


if __name__ == "__main__":
    main()
