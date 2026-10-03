from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

from nabu_protein.higher_order import (
    fnv1a32,
    model_diagnostics,
    parse_mutations,
)
from nabu_protein.metrics import evaluate_v83
from nabu_protein.v83 import NabuV83Model


def main(csv_path: str, frozen_result_path: str, out_path: str) -> None:
    data = pd.read_csv(csv_path)
    mutation_column = next(
        (c for c in ["mutant", "mutation", "mutations", "variant"] if c in data.columns),
        None,
    )
    truth_column = next(
        (c for c in ["DMS_score", "score", "fitness", "mean", "Fitness"] if c in data.columns),
        None,
    )
    if mutation_column is None or truth_column is None:
        raise RuntimeError(f"Could not identify required columns: {list(data.columns)}")

    data = data[[mutation_column, truth_column]].dropna().copy()
    data["mutation_set"] = data[mutation_column].map(parse_mutations)
    data = data[
        data["mutation_set"].map(
            lambda ms: 3 <= len(ms) <= 5 and "*" not in "".join(ms)
        )
    ].copy()
    data["candidate_id"] = data["mutation_set"].map(lambda ms: ":".join(ms))
    if data["candidate_id"].duplicated().any():
        raise RuntimeError("Duplicate canonical candidate IDs in external parity dataset.")

    data["bucket"] = data["candidate_id"].map(lambda cid: fnv1a32(cid) % 10)
    visible = data[data["bucket"] <= 6].copy()
    hidden = data[data["bucket"] >= 7].copy()

    package = NabuV83Model().fit(
        visible["mutation_set"].tolist(),
        visible[truth_column].to_numpy(dtype=float),
        visible["candidate_id"].astype(str).tolist(),
    )

    scored = package.score_candidates(
        hidden["mutation_set"].tolist(),
        hidden["candidate_id"].astype(str).tolist(),
    )
    eligible = hidden.merge(
        scored.drop(columns=["mutation_set"]),
        on="candidate_id",
        how="inner",
        validate="one_to_one",
    )
    eligible = eligible[eligible["scoreable"]].copy()

    frozen = json.loads(Path(frozen_result_path).read_text(encoding="utf-8"))

    package_metrics = {
        arm: evaluate_v83(eligible, arm, truth_column)
        for arm in [
            "B3_RAW_PAIR",
            "B4_CROSSFIT_TRIPLET",
            "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER",
            "V8_3_ADAPTIVE_ROUTER",
        ]
    }

    failures = []
    if package.router_decision["mode"] != frozen["oof"]["router_mode"]:
        failures.append(
            f"router mode {package.router_decision['mode']} != {frozen['oof']['router_mode']}"
        )

    for arm, actual_metrics in package_metrics.items():
        expected_metrics = frozen["metrics"][arm]
        for key, expected in expected_metrics.items():
            actual = actual_metrics[key]
            if isinstance(expected, int):
                if int(actual) != int(expected):
                    failures.append(f"{arm}.{key}: {actual} != {expected}")
            elif not np.isclose(
                float(actual),
                float(expected),
                rtol=1e-11,
                atol=1e-12,
            ):
                failures.append(f"{arm}.{key}: {actual} != {expected}")

    actual_top50 = eligible.sort_values(
        ["V8_3_ADAPTIVE_ROUTER", "candidate_id"],
        ascending=[False, True],
    ).head(50)["candidate_id"].astype(str).tolist()

    summary = {
        "dataset": Path(csv_path).name,
        "rows": int(len(data)),
        "visible": int(len(visible)),
        "hidden": int(len(hidden)),
        "eligible_hidden": int(len(eligible)),
        "router_mode": package.router_decision["mode"],
        "frozen_router_mode": frozen["oof"]["router_mode"],
        "model_diagnostics": model_diagnostics(package.model),
        "metrics": package_metrics,
        "top50_candidate_ids": actual_top50,
        "parity_passed": not failures,
        "failures": failures,
    }

    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(json.dumps(summary, indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("csv")
    parser.add_argument("frozen_result")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    main(args.csv, args.frozen_result, args.out)
