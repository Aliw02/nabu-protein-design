"""
Convenience benchmark CLI for the frozen NABU V8.3 core.

This command is deterministic but is NOT a substitute for a preregistered
Stage-A/Stage-B sealed scientific validation because the input CSV contains
both visible and hidden truth in one file.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from .higher_order import FROZEN_CORE_COMMIT, fnv1a32, parse_mutations
from .metrics import evaluate_v83
from .v83 import NabuV83Model


MUTATION_COLUMNS = ["mutant", "mutation", "mutations", "variant"]
FITNESS_COLUMNS = ["DMS_score", "score", "fitness", "mean", "Fitness"]


def find_column(columns, candidates):
    return next((name for name in candidates if name in columns), None)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Run a deterministic convenience benchmark with the frozen "
            "NABU V8.3 B2/B3/B4/B5 dual-objective router."
        )
    )
    parser.add_argument("csv", help="DMS CSV containing mutation and fitness columns")
    parser.add_argument(
        "--out",
        default="results_nabu_v83",
        help="Output directory (default: results_nabu_v83)",
    )
    parser.add_argument("--min-mutations", type=int, default=3)
    parser.add_argument("--max-mutations", type=int, default=5)
    args = parser.parse_args()

    data = pd.read_csv(args.csv)
    mutation_column = find_column(data.columns, MUTATION_COLUMNS)
    truth_column = find_column(data.columns, FITNESS_COLUMNS)
    if mutation_column is None or truth_column is None:
        raise SystemExit(
            "Could not identify mutation/fitness columns. "
            f"Observed columns: {list(data.columns)}"
        )

    data = data[[mutation_column, truth_column]].dropna().copy()
    data["mutation_set"] = data[mutation_column].map(parse_mutations)
    data = data[
        data["mutation_set"].map(
            lambda ms: (
                args.min_mutations <= len(ms) <= args.max_mutations
                and "*" not in "".join(ms)
            )
        )
    ].copy()

    data["candidate_id"] = data["mutation_set"].map(lambda ms: ":".join(ms))
    if data["candidate_id"].duplicated().any():
        count = int(data["candidate_id"].duplicated(keep=False).sum())
        raise SystemExit(
            f"Duplicate canonical candidate rows detected ({count}). "
            "Aggregate or deduplicate replicates before benchmarking."
        )

    data["bucket"] = data["candidate_id"].map(lambda cid: fnv1a32(cid) % 10)
    visible = data[data["bucket"] <= 6].copy()
    hidden = data[data["bucket"] >= 7].copy()

    if visible.empty or hidden.empty:
        raise SystemExit(
            f"Deterministic split is empty: visible={len(visible)} hidden={len(hidden)}"
        )

    model = NabuV83Model().fit(
        visible["mutation_set"].tolist(),
        visible[truth_column].to_numpy(dtype=float),
        visible["candidate_id"].astype(str).tolist(),
    )

    scored = model.score_candidates(
        hidden["mutation_set"].tolist(),
        hidden["candidate_id"].astype(str).tolist(),
    )
    hidden_scored = hidden.merge(
        scored.drop(columns=["mutation_set"]),
        on="candidate_id",
        how="left",
        validate="one_to_one",
    )
    eligible = hidden_scored[hidden_scored["scoreable"]].copy()
    if eligible.empty:
        raise SystemExit("No hidden candidates satisfy the frozen scoreability rule.")

    metrics = evaluate_v83(
        eligible,
        "V8_3_ADAPTIVE_ROUTER",
        truth_column,
    )

    ranked = eligible.sort_values(
        ["V8_3_ADAPTIVE_ROUTER", "candidate_id"],
        ascending=[False, True],
    )

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    ranked.drop(columns=["mutation_set"]).to_csv(
        out / "ranked_hidden_candidates.csv",
        index=False,
    )

    result = {
        "version": "NABU_V8_3_PACKAGE_BENCHMARK_V1",
        "frozen_core_commit": FROZEN_CORE_COMMIT,
        "scientific_sealed_validation": False,
        "dataset": str(args.csv),
        "mutation_column": mutation_column,
        "truth_column": truth_column,
        "total_filtered_rows": int(len(data)),
        "visible_rows": int(len(visible)),
        "hidden_rows": int(len(hidden)),
        "eligible_hidden_rows": int(len(eligible)),
        "router": model.router_decision,
        "metrics": metrics,
    }
    (out / "results.json").write_text(
        json.dumps(result, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
