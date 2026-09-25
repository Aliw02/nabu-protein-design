from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from nabu_protein.higher_order import fnv1a32, model_diagnostics, parse_mutations
from nabu_protein.v83 import NabuV83Model


def main(sealed_dir: str, expected_csv: str, expected_manifest: str, out_path: str) -> None:
    sealed = Path(sealed_dir)
    training_ids = pd.read_csv(sealed / "TRAINING_IDENTITIES.csv")
    training_labels = pd.read_csv(sealed / ".TRAINING_LABELS.sealed.csv")
    hidden = pd.read_csv(sealed / "HIDDEN_4_5_IDS.csv")

    manifest = json.loads(Path(expected_manifest).read_text(encoding="utf-8"))
    expected = pd.read_csv(expected_csv)

    ordered = training_ids.copy()
    ordered["_budget_hash"] = ordered["candidate_id"].map(
        lambda cid: fnv1a32("BUDGET|" + str(cid))
    )
    ordered = ordered.sort_values(
        ["_budget_hash", "candidate_id"],
        ascending=[True, True],
    ).reset_index(drop=True)

    count = int(manifest["budget_counts"]["40"])
    selected = ordered.head(count).copy()
    visible = selected.merge(
        training_labels[["candidate_id", "DMS_score"]],
        on="candidate_id",
        how="inner",
        validate="one_to_one",
    )
    if len(visible) != count:
        raise RuntimeError(f"Passive-40 label reveal mismatch: {len(visible)} != {count}")

    model = NabuV83Model().fit(
        visible["mutant"].map(parse_mutations).tolist(),
        visible["DMS_score"].to_numpy(dtype=float),
        visible["candidate_id"].astype(str).tolist(),
    )

    expected_mode = manifest["router_by_budget"]["passive"]["40"]["mode"]
    failures = []
    if model.router_decision["mode"] != expected_mode:
        failures.append(
            f"router mode {model.router_decision['mode']} != {expected_mode}"
        )

    actual = model.score_candidates(
        hidden["mutant"].map(parse_mutations).tolist(),
        hidden["candidate_id"].astype(str).tolist(),
    )
    actual = actual[actual["scoreable"]].copy()

    columns = {
        "P40_B3": "B3_RAW_PAIR",
        "P40_B5": "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER",
        "P40_V83": "V8_3_ADAPTIVE_ROUTER",
    }
    merged = expected[["candidate_id", *columns.keys()]].merge(
        actual[["candidate_id", *columns.values()]],
        on="candidate_id",
        how="inner",
        validate="one_to_one",
    )

    if len(merged) != len(expected) or len(actual) != len(expected):
        failures.append(
            f"eligible count mismatch actual={len(actual)} expected={len(expected)} merged={len(merged)}"
        )

    for expected_column, actual_column in columns.items():
        lhs = merged[actual_column].to_numpy(dtype=float)
        rhs = merged[expected_column].to_numpy(dtype=float)
        if not np.allclose(lhs, rhs, rtol=1e-13, atol=1e-13):
            diff = float(np.max(np.abs(lhs - rhs)))
            failures.append(f"{actual_column} max_abs_diff={diff}")

    actual_top50 = actual.sort_values(
        ["V8_3_ADAPTIVE_ROUTER", "candidate_id"],
        ascending=[False, True],
    ).head(50)["candidate_id"].astype(str).tolist()
    expected_top50 = expected.sort_values(
        ["P40_V83", "candidate_id"],
        ascending=[False, True],
    ).head(50)["candidate_id"].astype(str).tolist()
    if actual_top50 != expected_top50:
        failures.append("Passive-40 Top-50 candidate ordering differs.")

    result = {
        "version": "NABU_CR9114_PASSIVE40_PACKAGE_PARITY_V1",
        "visible_rows": int(len(visible)),
        "eligible_hidden_rows": int(len(actual)),
        "router_mode": model.router_decision["mode"],
        "expected_router_mode": expected_mode,
        "model_diagnostics": model_diagnostics(model.model),
        "top50_candidate_ids": actual_top50,
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
    parser.add_argument("sealed_dir")
    parser.add_argument("expected_csv")
    parser.add_argument("expected_manifest")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    main(
        args.sealed_dir,
        args.expected_csv,
        args.expected_manifest,
        args.out,
    )
