from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from nabu_protein.higher_order import parse_mutations
from nabu_protein.v83 import NabuV83Model


ROOT = Path(__file__).resolve().parents[1]


CASES = [
    {
        "name": "eqFP611",
        "visible": "eqfp611_sealed_input/VISIBLE.csv",
        "hidden": "eqfp611_sealed_input/HIDDEN_IDS.csv",
        "expected": "nabu_v8_3_eqfp611_stage_a/ALL_HIDDEN_SCORES.csv",
        "manifest": "nabu_v8_3_eqfp611_stage_a/STAGE_A_MANIFEST.json",
        "expected_mode": "GLOBAL_HIGHER_ORDER",
    },
    {
        "name": "CreiLOV",
        "visible": "creilov_sealed_input/VISIBLE.csv",
        "hidden": "creilov_sealed_input/HIDDEN_IDS.csv",
        "expected": "nabu_v8_3_creilov_stage_a/ALL_HIDDEN_SCORES.csv",
        "manifest": "nabu_v8_3_creilov_stage_a/STAGE_A_MANIFEST.json",
        "expected_mode": "B3_PROTECTED_NO_HIGHER_ORDER",
    },
]


@pytest.mark.parametrize("case", CASES, ids=[case["name"] for case in CASES])
def test_packaged_v83_reproduces_frozen_stage_a_scores(case):
    visible = pd.read_csv(ROOT / case["visible"])
    hidden = pd.read_csv(ROOT / case["hidden"])
    expected = pd.read_csv(ROOT / case["expected"])
    manifest = json.loads((ROOT / case["manifest"]).read_text(encoding="utf-8"))

    visible_sets = visible["mutant"].map(parse_mutations).tolist()
    hidden_sets = hidden["mutant"].map(parse_mutations).tolist()

    model = NabuV83Model().fit(
        visible_sets,
        visible["DMS_score"].to_numpy(dtype=float),
        visible["candidate_id"].astype(str).tolist(),
    )

    assert model.router_decision["mode"] == case["expected_mode"]
    assert manifest["real_router"]["mode"] == case["expected_mode"]

    actual = model.score_candidates(
        hidden_sets,
        hidden["candidate_id"].astype(str).tolist(),
    )
    actual = actual[actual["scoreable"]].copy()

    compare_columns = [
        "B3_RAW_PAIR",
        "B4_CROSSFIT_TRIPLET",
        "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER",
        "V8_3_ADAPTIVE_ROUTER",
        "triplet_delta",
        "quartet_delta",
        "triplet_confidence",
        "quartet_confidence",
    ]

    merged = expected[
        ["candidate_id", *compare_columns]
    ].merge(
        actual[["candidate_id", *compare_columns]],
        on="candidate_id",
        suffixes=("_expected", "_actual"),
        how="inner",
        validate="one_to_one",
    )

    assert len(merged) == len(expected) == len(actual)

    for column in compare_columns:
        np.testing.assert_allclose(
            merged[f"{column}_actual"].to_numpy(dtype=float),
            merged[f"{column}_expected"].to_numpy(dtype=float),
            rtol=1e-13,
            atol=1e-13,
        )

    expected_top50 = expected.sort_values(
        ["V8_3_ADAPTIVE_ROUTER", "candidate_id"],
        ascending=[False, True],
    ).head(50)["candidate_id"].astype(str).tolist()

    actual_top50 = actual.sort_values(
        ["V8_3_ADAPTIVE_ROUTER", "candidate_id"],
        ascending=[False, True],
    ).head(50)["candidate_id"].astype(str).tolist()

    assert actual_top50 == expected_top50
