from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import pandas as pd

from campaign import (
    CampaignSimulator,
    IdentityHashPolicy,
    VirtualAssayOracle,
    identity_pool_from_frame,
    truth_frame_from_frame,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run the Phase-2A.0 leak-safe campaign simulator on the RhlA "
            "training pool. This is an engineering smoke test, not a new "
            "scientific benchmark."
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
            "2A0_campaign_simulator/results/rhla_smoke"
        ),
    )
    parser.add_argument("--seed-fraction", type=float, default=0.05)
    parser.add_argument("--rounds", type=int, default=3)
    parser.add_argument("--batch-size", type=int, default=16)
    args = parser.parse_args()

    if not (0.0 < args.seed_fraction < 1.0):
        raise SystemExit("--seed-fraction must be between 0 and 1.")
    if args.rounds < 1:
        raise SystemExit("--rounds must be positive.")
    if args.batch_size < 1:
        raise SystemExit("--batch-size must be positive.")

    source = pd.read_csv(args.input)
    identity_pool = identity_pool_from_frame(source)
    truth = truth_frame_from_frame(source)
    oracle = VirtualAssayOracle(truth)

    out = Path(args.out)
    simulator = CampaignSimulator(
        identity_pool=identity_pool,
        oracle=oracle,
        output_dir=out,
    )

    seed_count = int(math.ceil(len(identity_pool) * args.seed_fraction))
    simulator.bootstrap(seed_count=seed_count, min_per_fold=1)

    policy = IdentityHashPolicy()
    for _ in range(args.rounds):
        if simulator.measurements_spent >= len(identity_pool):
            break
        simulator.run_round(policy, batch_size=args.batch_size)

    summary = simulator.write_campaign_manifest()
    summary["engineering_smoke_only"] = True
    summary["scientific_claim"] = False
    summary["input"] = str(args.input)
    summary["seed_fraction_requested"] = float(args.seed_fraction)
    summary["batch_size"] = int(args.batch_size)

    (out / "SMOKE_SUMMARY.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
