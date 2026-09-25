from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
PHASE2 = HERE.parent
PHASE2A0 = PHASE2 / "2A_closed_loop_acquisition" / "2A0_campaign_simulator"
PHASE2A1 = PHASE2 / "2A_closed_loop_acquisition" / "2A1_baselines"
PHASE2B = PHASE2 / "2B_candidate_assembly"
for path in (PHASE2A0, PHASE2A1, PHASE2B, HERE):
    sys.path.insert(0, str(path))

from campaign import identity_pool_from_frame  # noqa: E402
from nabu_protein.higher_order import fnv1a32  # noqa: E402
from gb1_adapter import adapt_gb1_frame, gb1_vocabulary  # noqa: E402
from combined_loop import (  # noqa: E402
    DesignCampaign,
    EvidenceMaturityGate,
    FullLandscapeEvaluator,
    FullLandscapeOracle,
    make_policy,
    selection_sha256,
)


RESERVOIR_NAMESPACE = "NABU_PHASE2C_GB1_RESERVOIR_V1"


def deterministic_reservoir(
    full_frame: pd.DataFrame,
    reservoir_size: int,
) -> pd.DataFrame:
    if len(full_frame) < reservoir_size:
        raise ValueError("GB1 universe smaller than requested reservoir.")

    work = full_frame[
        ["candidate_id", "mutant", "mutation_set"]
    ].copy()
    work["_hash"] = work["candidate_id"].astype(str).map(
        lambda candidate_id: fnv1a32(
            f"{RESERVOIR_NAMESPACE}|{candidate_id}"
        )
    )
    work = work.sort_values(
        ["_hash", "candidate_id"],
        ascending=[True, True],
    ).head(reservoir_size)
    work = work.drop(columns=["_hash"]).reset_index(drop=True)

    parser_frame = work[["candidate_id", "mutant"]].copy()
    parsed = identity_pool_from_frame(parser_frame)
    parsed["mutation_set"] = work["mutation_set"].tolist()
    return parsed


def json_default(value):
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    raise TypeError(f"Not JSON serializable: {type(value)!r}")


def next_decision(
    campaign: DesignCampaign,
    gate: EvidenceMaturityGate,
    pre_gate_policy,
    batch_size: int,
    use_assembly: bool,
):
    diagnostics, mature = gate.probe(
        model=campaign.model,
        measured_mutation_sets=campaign.measured["mutation_set"].tolist(),
        measured_ids=campaign.measured["candidate_id"].astype(str).tolist(),
    )

    if use_assembly and diagnostics.ready:
        selected = (
            mature.head(batch_size)["candidate_id"]
            .astype(str)
            .tolist()
        )
        return selected, "assembly", diagnostics

    selected = campaign.acquisition_selection(
        pre_gate_policy,
        batch_size=batch_size,
    )
    return selected, "acquisition", diagnostics


def hidden_truth_perturbation_check(
    full_frame: pd.DataFrame,
    reservoir: pd.DataFrame,
    seed_count: int,
    gate: EvidenceMaturityGate,
    policy_name: str,
    batch_size: int,
    use_assembly: bool,
) -> bool:
    left = DesignCampaign(full_frame, reservoir)
    left.bootstrap(seed_count)

    perturbed = full_frame.copy()
    seed_ids = set(left.bootstrap_ids)
    mask = ~perturbed["candidate_id"].astype(str).isin(seed_ids)
    perturbed.loc[mask, "DMS_score"] = (
        -1000000.0
        - np.arange(int(mask.sum()), dtype=float)
    )

    right = DesignCampaign(perturbed, reservoir)
    right.bootstrap(seed_count)

    left_policy = make_policy(policy_name, full_frame)
    right_policy = make_policy(policy_name, full_frame)

    left_ids, left_source, left_gate = next_decision(
        left, gate, left_policy, batch_size, use_assembly
    )
    right_ids, right_source, right_gate = next_decision(
        right, gate, right_policy, batch_size, use_assembly
    )

    return bool(
        left_ids == right_ids
        and left_source == right_source
        and left_gate.ready == right_gate.ready
        and left_gate.mature_count == right_gate.mature_count
    )


def state_matched_ablation(
    campaign: DesignCampaign,
    gate: EvidenceMaturityGate,
    acquisition_policy,
    batch_size: int,
) -> dict | None:
    diagnostics, mature = gate.probe(
        model=campaign.model,
        measured_mutation_sets=campaign.measured["mutation_set"].tolist(),
        measured_ids=campaign.measured["candidate_id"].astype(str).tolist(),
    )
    if not diagnostics.ready:
        return None

    assembly_ids = (
        mature.head(batch_size)["candidate_id"].astype(str).tolist()
    )
    acquisition_ids = campaign.acquisition_selection(
        acquisition_policy,
        batch_size=batch_size,
    )

    frozen = {
        "measurements_spent": campaign.measurements_spent,
        "assembly_ids": assembly_ids,
        "assembly_sha256": selection_sha256(assembly_ids),
        "acquisition_ids": acquisition_ids,
        "acquisition_sha256": selection_sha256(acquisition_ids),
        "truth_exposed_before_both_freezes": False,
        "gate_mature_count": diagnostics.mature_count,
    }

    assembly_truth = campaign.oracle.reveal(assembly_ids)
    acquisition_truth = campaign.oracle.reveal(acquisition_ids)

    cutoff = float(
        np.quantile(
            campaign.full_frame["DMS_score"].to_numpy(dtype=float),
            0.99,
        )
    )
    frozen["post_freeze_evaluation"] = {
        "assembly_best": float(
            assembly_truth["assay_value"].max()
        ),
        "acquisition_best": float(
            acquisition_truth["assay_value"].max()
        ),
        "assembly_mean": float(
            assembly_truth["assay_value"].mean()
        ),
        "acquisition_mean": float(
            acquisition_truth["assay_value"].mean()
        ),
        "assembly_top1_hits": int(
            (assembly_truth["assay_value"] >= cutoff).sum()
        ),
        "acquisition_top1_hits": int(
            (acquisition_truth["assay_value"] >= cutoff).sum()
        ),
    }
    return frozen


def run_condition(
    condition_name: str,
    pre_gate_policy_name: str,
    use_assembly: bool,
    full_frame: pd.DataFrame,
    reservoir: pd.DataFrame,
    vocabulary,
    output_dir: Path,
    seed_count: int,
    target_count: int,
    batch_size: int,
    probe_count: int,
) -> dict:
    condition_dir = output_dir / condition_name
    campaign = DesignCampaign(
        full_frame=full_frame,
        reservoir=reservoir,
        output_dir=condition_dir / "campaign",
    )
    bootstrap = campaign.bootstrap(seed_count)

    policy = make_policy(pre_gate_policy_name, full_frame)
    gate = EvidenceMaturityGate(
        vocabulary=vocabulary,
        reservoir_ids=reservoir["candidate_id"].astype(str),
        assayable_ids=full_frame["candidate_id"].astype(str),
        batch_size=batch_size,
        min_main_support=2,
        min_pair_support=2,
        beam_width=256,
        probe_count=int(probe_count),
        target_order=3,
    )

    evaluator = FullLandscapeEvaluator(
        full_frame=full_frame,
        bootstrap_count=campaign.measurements_spent,
    )
    curve = [evaluator.checkpoint(campaign)]
    gate_rows = []
    first_gate_open_measurement = None
    matched_ablation = None

    while campaign.measurements_spent < target_count:
        remaining = target_count - campaign.measurements_spent
        actual_batch = min(batch_size, remaining)

        selected, source, diagnostics = next_decision(
            campaign=campaign,
            gate=gate,
            pre_gate_policy=policy,
            batch_size=actual_batch,
            use_assembly=use_assembly,
        )

        gate_rows.append(
            {
                "measurements_spent": campaign.measurements_spent,
                "ready": diagnostics.ready,
                "mature_count": diagnostics.mature_count,
                "probe_count": diagnostics.probe_count,
                "selected_source": source,
            }
        )

        if diagnostics.ready and first_gate_open_measurement is None:
            first_gate_open_measurement = campaign.measurements_spent
            if condition_name == "acquisition_then_gated_assembly":
                matched_ablation = state_matched_ablation(
                    campaign=campaign,
                    gate=gate,
                    acquisition_policy=make_policy(
                        "historical_50_50",
                        full_frame,
                    ),
                    batch_size=actual_batch,
                )

        campaign.reveal_frozen_batch(
            selected,
            source=source,
            gate=diagnostics,
        )
        curve.append(evaluator.checkpoint(campaign))

    pd.DataFrame(curve).to_csv(
        condition_dir / "DISCOVERY_CURVE.csv",
        index=False,
    )
    pd.DataFrame(gate_rows).to_csv(
        condition_dir / "GATE_TRAJECTORY.csv",
        index=False,
    )
    if matched_ablation is not None:
        (condition_dir / "STATE_MATCHED_ABLATION.json").write_text(
            json.dumps(
                matched_ablation,
                indent=2,
                default=json_default,
            ),
            encoding="utf-8",
        )

    summary = evaluator.summarize(curve)
    summary.update(
        {
            "version": "NABU_PHASE2C_CONDITION_V1",
            "condition": condition_name,
            "pre_gate_policy": pre_gate_policy_name,
            "assembly_enabled": bool(use_assembly),
            "seed_count": int(seed_count),
            "target_count": int(target_count),
            "batch_size": int(batch_size),
            "bootstrap_sha256": bootstrap[
                "selection_sha256_before_label_reveal"
            ],
            "first_gate_open_measurement": first_gate_open_measurement,
            "assembly_rounds": int(
                sum(
                    row["selected_source"] == "assembly"
                    for row in gate_rows
                )
            ),
            "gate_ready_rounds": int(
                sum(bool(row["ready"]) for row in gate_rows)
            ),
        }
    )

    (condition_dir / "CONDITION_SUMMARY.json").write_text(
        json.dumps(summary, indent=2, default=json_default),
        encoding="utf-8",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run preregistered NABU Phase-2C GB1 POC."
    )
    parser.add_argument("source_csv")
    parser.add_argument(
        "--out",
        default=(
            "phase2/2C_assembly_plus_acquisition/"
            "results/gb1_gated_v1"
        ),
    )
    parser.add_argument("--reservoir-size", type=int, default=2048)
    parser.add_argument("--seed-count", type=int, default=64)
    parser.add_argument("--target-count", type=int, default=256)
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--probe-count", type=int, default=64)
    args = parser.parse_args()

    raw = pd.read_csv(args.source_csv)
    full_frame = adapt_gb1_frame(raw)
    reservoir = deterministic_reservoir(
        full_frame,
        reservoir_size=args.reservoir_size,
    )
    vocabulary = gb1_vocabulary()

    if not (
        0 < args.seed_count < args.target_count <= args.reservoir_size
    ):
        raise SystemExit(
            "Require 0 < seed-count < target-count <= reservoir-size."
        )

    gate = EvidenceMaturityGate(
        vocabulary=vocabulary,
        reservoir_ids=reservoir["candidate_id"].astype(str),
        assayable_ids=full_frame["candidate_id"].astype(str),
        batch_size=args.batch_size,
        min_main_support=2,
        min_pair_support=2,
        beam_width=256,
        probe_count=int(args.probe_count),
        target_order=3,
    )

    leakage_checks = {
        "random_only": hidden_truth_perturbation_check(
            full_frame,
            reservoir,
            args.seed_count,
            gate,
            "random",
            args.batch_size,
            False,
        ),
        "acquisition_only": hidden_truth_perturbation_check(
            full_frame,
            reservoir,
            args.seed_count,
            gate,
            "historical_50_50",
            args.batch_size,
            False,
        ),
        "random_then_gated_assembly": hidden_truth_perturbation_check(
            full_frame,
            reservoir,
            args.seed_count,
            gate,
            "random",
            args.batch_size,
            True,
        ),
        "acquisition_then_gated_assembly": hidden_truth_perturbation_check(
            full_frame,
            reservoir,
            args.seed_count,
            gate,
            "historical_50_50",
            args.batch_size,
            True,
        ),
    }
    if not all(leakage_checks.values()):
        raise RuntimeError(
            f"Hidden-truth perturbation failed: {leakage_checks}"
        )

    matrix = [
        ("random_only", "random", False),
        ("acquisition_only", "historical_50_50", False),
        ("random_then_gated_assembly", "random", True),
        (
            "acquisition_then_gated_assembly",
            "historical_50_50",
            True,
        ),
    ]

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    summaries = []
    for condition_name, policy_name, use_assembly in matrix:
        print(f"[2C] START {condition_name}")
        summary = run_condition(
            condition_name=condition_name,
            pre_gate_policy_name=policy_name,
            use_assembly=use_assembly,
            full_frame=full_frame,
            reservoir=reservoir,
            vocabulary=vocabulary,
            output_dir=out,
            seed_count=args.seed_count,
            target_count=args.target_count,
            batch_size=args.batch_size,
            probe_count=args.probe_count,
        )
        summaries.append(summary)
        print(
            f"[2C] DONE {condition_name} "
            f"best={summary['final']['best_true_fitness']:.6f} "
            f"top1={summary['final']['cumulative_top1_hits']} "
            f"auc={summary['post_bootstrap_normalized_discovery_auc']:.6f} "
            f"assembly_rounds={summary['assembly_rounds']}"
        )

    bootstrap_hashes = {
        summary["bootstrap_sha256"]
        for summary in summaries
    }
    if len(bootstrap_hashes) != 1:
        raise RuntimeError(
            "2C conditions did not share the exact same bootstrap."
        )

    rows = []
    for summary in summaries:
        final = summary["final"]
        rows.append(
            {
                "condition": summary["condition"],
                "post_bootstrap_normalized_discovery_auc": summary[
                    "post_bootstrap_normalized_discovery_auc"
                ],
                "final_best_true_fitness": final["best_true_fitness"],
                "final_regret": final["best_so_far_regret"],
                "final_top1_hits": final["cumulative_top1_hits"],
                "assembly_top1_hits": final["assembly_top1_hits"],
                "first_top1_measurement": summary[
                    "first_top1_measurement"
                ],
                "first_gate_open_measurement": summary[
                    "first_gate_open_measurement"
                ],
                "assembly_rounds": summary["assembly_rounds"],
            }
        )
    comparison = pd.DataFrame(rows)
    comparison.to_csv(out / "PHASE2C_COMPARISON.csv", index=False)

    by_name = comparison.set_index("condition")
    acquisition = by_name.loc["acquisition_only"]
    combined = by_name.loc["acquisition_then_gated_assembly"]

    scientific_checks = {
        "auc_strictly_improves": bool(
            combined["post_bootstrap_normalized_discovery_auc"]
            > acquisition["post_bootstrap_normalized_discovery_auc"]
        ),
        "final_regret_nonworse": bool(
            combined["final_regret"] <= acquisition["final_regret"]
        ),
        "top1_hits_nonworse": bool(
            combined["final_top1_hits"] >= acquisition["final_top1_hits"]
        ),
        "assembly_activated": bool(
            combined["assembly_rounds"] >= 1
        ),
    }
    scientific_pass = all(scientific_checks.values())

    manifest = {
        "version": (
            "NABU_PHASE2C_GB1_GATED_V2"
            if int(args.probe_count) == 256
            else "NABU_PHASE2C_GB1_GATED_V1"
        ),
        "scientific_claim": False,
        "stage": "2C_ASSEMBLY_PLUS_ACQUISITION_DEVELOPMENT",
        "development_dataset": "GB1_Wu2016",
        "gb1_reference_sequence": (
            "MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE"
        ),
        "mutable_positions": [39, 40, 41, 54],
        "wt_genotype": "VDGV",
        "assay_universe_size": int(len(full_frame)),
        "reservoir_size": int(len(reservoir)),
        "reservoir_namespace": RESERVOIR_NAMESPACE,
        "seed_count": int(args.seed_count),
        "target_count": int(args.target_count),
        "batch_size": int(args.batch_size),
        "gate": {
            "target_order": 3,
            "probe_count": int(args.probe_count),
            "beam_width": 256,
            "min_main_support": 2,
            "min_pair_support": 2,
            "requires_outside_reservoir": True,
            "requires_identity_in_assay_universe": True,
        },
        "hidden_truth_perturbation_checks": leakage_checks,
        "comparison": rows,
        "scientific_development_checks": scientific_checks,
        "scientific_development_pass": scientific_pass,
        "interpretation_rule": (
            "This GB1 run is development evidence only. "
            "Do not retune V1 against these outcomes and call the rerun blind. "
            "Final Phase-2 PASS remains reserved for 2D."
        ),
    }
    (out / "PHASE2C_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, default=json_default),
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2, default=json_default))


if __name__ == "__main__":
    main()
