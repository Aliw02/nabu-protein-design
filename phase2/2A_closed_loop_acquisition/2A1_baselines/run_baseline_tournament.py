from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
PHASE2A0 = HERE.parent / "2A0_campaign_simulator"
sys.path.insert(0, str(PHASE2A0))

from campaign import (  # noqa: E402
    CampaignSimulator,
    LABEL_COLUMN,
    VirtualAssayOracle,
    identity_pool_from_frame,
    truth_frame_from_frame,
)
from baselines import build_baselines  # noqa: E402


def json_default(value):
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    raise TypeError(f"Not JSON serializable: {type(value)!r}")


def pairwise_jaccard_distance(mutation_sets: list[tuple[str, ...]]) -> float | None:
    if len(mutation_sets) < 2:
        return None

    distances = []
    for i in range(len(mutation_sets)):
        left = set(mutation_sets[i])
        for j in range(i + 1, len(mutation_sets)):
            right = set(mutation_sets[j])
            union = left | right
            similarity = len(left & right) / len(union) if union else 1.0
            distances.append(1.0 - similarity)

    return float(np.mean(distances)) if distances else None


class CampaignEvaluator:
    """Full-truth evaluator. Never passed to an acquisition policy."""

    def __init__(
        self,
        identity_pool: pd.DataFrame,
        truth: pd.DataFrame,
        bootstrap_count: int,
    ):
        joined = identity_pool[
            ["candidate_id", "mutation_set"]
        ].merge(
            truth,
            on="candidate_id",
            how="inner",
            validate="one_to_one",
        )
        if len(joined) != len(identity_pool):
            raise RuntimeError("Evaluator truth does not cover identity pool.")

        self.truth_by_id = dict(
            zip(
                joined["candidate_id"].astype(str),
                joined[LABEL_COLUMN].astype(float),
            )
        )
        self.mutation_by_id = dict(
            zip(
                joined["candidate_id"].astype(str),
                joined["mutation_set"],
            )
        )
        values = joined[LABEL_COLUMN].to_numpy(dtype=float)
        self.global_best = float(np.max(values))
        self.global_min = float(np.min(values))
        self.top1_cutoff = float(np.quantile(values, 0.99))
        self.bootstrap_count = int(bootstrap_count)

    def checkpoint(self, simulator: CampaignSimulator) -> dict:
        measured_ids = simulator.measured[
            "candidate_id"
        ].astype(str).tolist()
        measured_values = np.array(
            [self.truth_by_id[candidate_id] for candidate_id in measured_ids],
            dtype=float,
        )
        best = float(np.max(measured_values))
        top_ids = [
            candidate_id
            for candidate_id in measured_ids
            if self.truth_by_id[candidate_id] >= self.top1_cutoff
        ]

        if simulator.model is None:
            scoreable_coverage = None
            router_mode = None
            oof_spearman = None
        else:
            view = simulator.policy_view()
            if len(view.candidates):
                scoreable_coverage = float(
                    view.candidates["scoreable"].astype(bool).mean()
                )
            else:
                scoreable_coverage = 1.0
            router_mode = str(simulator.model.router_decision["mode"])
            oof_spearman = float(
                simulator.model.router_decision["B3_oof_spearman"]
            )

        denominator = self.global_best - self.global_min
        normalized_best = (
            1.0
            if denominator <= 0
            else float((best - self.global_min) / denominator)
        )

        return {
            "measurements_spent": int(simulator.measurements_spent),
            "best_true_fitness": best,
            "normalized_best_fitness": normalized_best,
            "best_so_far_regret": float(self.global_best - best),
            "cumulative_top1_hits": int(len(top_ids)),
            "top1_mean_pairwise_jaccard_distance": pairwise_jaccard_distance(
                [self.mutation_by_id[candidate_id] for candidate_id in top_ids]
            ),
            "scoreable_coverage_unmeasured": scoreable_coverage,
            "router_mode": router_mode,
            "visible_b3_oof_spearman": oof_spearman,
        }

    def summarize(
        self,
        curve: list[dict],
        bootstrap_ids: set[str],
        final_measured_ids: list[str],
    ) -> dict:
        frame = pd.DataFrame(curve).sort_values("measurements_spent")
        post = frame[
            frame["measurements_spent"] >= self.bootstrap_count
        ].copy()

        if len(post) < 2:
            auc = float(post["normalized_best_fitness"].iloc[0])
        else:
            x = post["measurements_spent"].to_numpy(dtype=float)
            y = post["normalized_best_fitness"].to_numpy(dtype=float)
            x_span = float(x[-1] - x[0])
            auc = (
                float(y[-1])
                if x_span <= 0
                else float(np.trapezoid(y, x) / x_span)
            )

        top1_rows = frame[frame["cumulative_top1_hits"] > 0]
        first_top1 = (
            None
            if top1_rows.empty
            else int(top1_rows.iloc[0]["measurements_spent"])
        )

        acquired_ids = [
            candidate_id
            for candidate_id in final_measured_ids
            if candidate_id not in bootstrap_ids
        ]
        acquired_top1 = sum(
            self.truth_by_id[candidate_id] >= self.top1_cutoff
            for candidate_id in acquired_ids
        )

        bootstrap_best = max(
            self.truth_by_id[candidate_id]
            for candidate_id in bootstrap_ids
        )
        acquired_best = (
            None
            if not acquired_ids
            else max(self.truth_by_id[candidate_id] for candidate_id in acquired_ids)
        )

        final = frame.iloc[-1].to_dict()
        return {
            "global_best_true_fitness": self.global_best,
            "top1_cutoff": self.top1_cutoff,
            "bootstrap_best_true_fitness": float(bootstrap_best),
            "best_acquired_after_bootstrap": (
                None if acquired_best is None else float(acquired_best)
            ),
            "acquired_top1_hits_after_bootstrap": int(acquired_top1),
            "first_top1_measurement": first_top1,
            "post_bootstrap_normalized_discovery_auc": auc,
            "final": final,
        }


def run_policy(
    policy,
    identity_pool: pd.DataFrame,
    truth: pd.DataFrame,
    output_dir: Path,
    seed_count: int,
    target_count: int,
    batch_size: int,
) -> dict:
    baseline_dir = output_dir / policy.name
    simulator = CampaignSimulator(
        identity_pool=identity_pool,
        oracle=VirtualAssayOracle(truth),
        output_dir=baseline_dir / "campaign",
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

    while simulator.measurements_spent < target_count:
        remaining = target_count - simulator.measurements_spent
        simulator.run_round(
            policy,
            batch_size=min(batch_size, remaining),
        )
        curve.append(evaluator.checkpoint(simulator))

    curve_frame = pd.DataFrame(curve)
    curve_frame.to_csv(
        baseline_dir / "DISCOVERY_CURVE.csv",
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
            "version": "NABU_PHASE2A1_BASELINE_RESULT_V1",
            "baseline": policy.name,
            "seed_count": int(seed_count),
            "target_count": int(target_count),
            "batch_size": int(batch_size),
            "rounds_completed": int(len(simulator.round_manifests)),
            "bootstrap_selection_sha256": bootstrap[
                "selection_sha256_before_label_reveal"
            ],
            "scientific_scope": (
                "RhlA training-pool virtual campaign baseline comparison; "
                "not an untouched-transfer claim."
            ),
        }
    )

    (baseline_dir / "BASELINE_SUMMARY.json").write_text(
        json.dumps(summary, indent=2, default=json_default),
        encoding="utf-8",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the preregistered Phase-2A.1 baseline tournament."
    )
    parser.add_argument(
        "--input",
        default="rhla_sample_efficiency_sealed_input/TRAINING_POOL.csv",
    )
    parser.add_argument(
        "--out",
        default=(
            "phase2/2A_closed_loop_acquisition/"
            "2A1_baselines/results/rhla_baseline_tournament"
        ),
    )
    parser.add_argument("--seed-fraction", type=float, default=0.05)
    parser.add_argument("--target-fraction", type=float, default=0.20)
    parser.add_argument("--batch-size", type=int, default=16)
    args = parser.parse_args()

    if not (0.0 < args.seed_fraction < args.target_fraction <= 1.0):
        raise SystemExit(
            "Require 0 < seed-fraction < target-fraction <= 1."
        )
    if args.batch_size < 1:
        raise SystemExit("--batch-size must be positive.")

    source = pd.read_csv(args.input)
    identity_pool = identity_pool_from_frame(source)
    truth = truth_frame_from_frame(source)

    seed_count = int(math.ceil(len(identity_pool) * args.seed_fraction))
    target_count = int(math.ceil(len(identity_pool) * args.target_fraction))

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    summaries = []
    for policy in build_baselines(identity_pool):
        print(f"[2A1] START {policy.name}")
        summary = run_policy(
            policy=policy,
            identity_pool=identity_pool,
            truth=truth,
            output_dir=out,
            seed_count=seed_count,
            target_count=target_count,
            batch_size=args.batch_size,
        )
        summaries.append(summary)
        print(
            f"[2A1] DONE {policy.name} "
            f"best={summary['final']['best_true_fitness']:.6f} "
            f"top1={summary['final']['cumulative_top1_hits']} "
            f"auc={summary['post_bootstrap_normalized_discovery_auc']:.6f}"
        )

    bootstrap_hashes = {
        summary["bootstrap_selection_sha256"]
        for summary in summaries
    }
    if len(bootstrap_hashes) != 1:
        raise RuntimeError("Baselines did not share the exact same bootstrap.")

    comparison_rows = []
    for summary in summaries:
        final = summary["final"]
        comparison_rows.append(
            {
                "baseline": summary["baseline"],
                "seed_count": summary["seed_count"],
                "target_count": summary["target_count"],
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
            }
        )

    comparison = pd.DataFrame(comparison_rows)
    comparison.to_csv(out / "BASELINE_COMPARISON.csv", index=False)

    manifest = {
        "version": "NABU_PHASE2A1_BASELINE_TOURNAMENT_V1",
        "scientific_claim": False,
        "stage": "2A1_BASELINE_CHARACTERIZATION",
        "input": str(args.input),
        "pool_size": int(len(identity_pool)),
        "seed_fraction": float(args.seed_fraction),
        "target_fraction": float(args.target_fraction),
        "seed_count": int(seed_count),
        "target_count": int(target_count),
        "batch_size": int(args.batch_size),
        "shared_bootstrap_selection_sha256": next(iter(bootstrap_hashes)),
        "baselines": [summary["baseline"] for summary in summaries],
        "comparison": comparison_rows,
        "interpretation_rule": (
            "Do not select a Phase-2 winner or tune a new controller from "
            "this single landscape alone. Use these results to establish "
            "controls for later multi-landscape work."
        ),
    }
    (out / "TOURNAMENT_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, default=json_default),
        encoding="utf-8",
    )

    print(json.dumps(manifest, indent=2, default=json_default))


if __name__ == "__main__":
    main()
