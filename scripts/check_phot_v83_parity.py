from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from nabu_protein.higher_order import fnv1a32, parse_mutations
from nabu_protein.metrics import evaluate_v83
from nabu_protein.v83 import NabuV83Model


def main(csv_path: str, frozen_result_path: str, out_path: str) -> None:
    data = pd.read_csv(csv_path)
    mutation_column = next(
        (c for c in ["mutant", "mutation", "mutations"] if c in data.columns),
        None,
    )
    truth_column = next(
        (c for c in ["DMS_score", "score", "fitness", "mean"] if c in data.columns),
        None,
    )
    if mutation_column is None or truth_column is None:
        raise RuntimeError(f"Required PHOT columns not found: {list(data.columns)}")

    data = data[[mutation_column, truth_column]].dropna().copy()
    data["mutation_set"] = data[mutation_column].map(parse_mutations)
    data = data[
        data["mutation_set"].map(
            lambda ms: 3 <= len(ms) <= 5 and "*" not in "".join(ms)
        )
    ].copy()
    data["candidate_id"] = data["mutation_set"].map(lambda ms: ":".join(ms))
    if data["candidate_id"].duplicated().any():
        raise RuntimeError("Duplicate canonical PHOT candidate IDs.")

    data["bucket"] = data["candidate_id"].map(lambda cid: fnv1a32(cid) % 10)
    visible = data[data["bucket"] <= 6].copy()
    hidden = data[data["bucket"] >= 7].copy()

    model = NabuV83Model().fit(
        visible["mutation_set"].tolist(),
        visible[truth_column].to_numpy(dtype=float),
        visible["candidate_id"].astype(str).tolist(),
    )
    scored = model.score_candidates(
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

    frozen_root = json.loads(Path(frozen_result_path).read_text(encoding="utf-8"))
    frozen = next(
        row
        for row in frozen_root["datasets"]
        if str(row["dataset"]).startswith("PHOT_CHLRE_Chen_2023")
    )

    failures = []
    if model.router_decision["mode"] != frozen["router_mode"]:
        failures.append(
            f"router mode {model.router_decision['mode']} != {frozen['router_mode']}"
        )

    actual_metrics = {
        arm: evaluate_v83(eligible, arm, truth_column)
        for arm in [
            "B3_RAW_PAIR",
            "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER",
            "V8_3_ADAPTIVE_ROUTER",
        ]
    }

    for arm, actual in actual_metrics.items():
        expected = frozen["metrics"][arm]
        for key, expected_value in expected.items():
            if key == "top50_top1pct_enrichment":
                actual_value = float(actual["top50_top1pct_hits"]) / 0.5
            else:
                actual_value = actual[key]
            if isinstance(expected_value, int):
                if int(actual_value) != int(expected_value):
                    failures.append(f"{arm}.{key}: {actual_value} != {expected_value}")
            elif not np.isclose(
                float(actual_value),
                float(expected_value),
                rtol=1e-11,
                atol=1e-12,
            ):
                failures.append(f"{arm}.{key}: {actual_value} != {expected_value}")

    result = {
        "version": "NABU_PHOT_PACKAGE_PARITY_V1",
        "rows": int(len(data)),
        "visible": int(len(visible)),
        "hidden": int(len(hidden)),
        "eligible_hidden": int(len(eligible)),
        "router_mode": model.router_decision["mode"],
        "expected_router_mode": frozen["router_mode"],
        "metrics": actual_metrics,
        "parity_passed": not failures,
        "failures": failures,
    }

    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))

    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("csv")
    parser.add_argument("frozen_result")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    main(args.csv, args.frozen_result, args.out)
