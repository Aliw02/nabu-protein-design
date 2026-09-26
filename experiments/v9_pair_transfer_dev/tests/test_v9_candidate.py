from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import numpy as np


RUNNER = (
    Path(__file__).resolve().parents[1]
    / "run_v9_candidate.py"
)
SPEC = importlib.util.spec_from_file_location("run_v9_candidate", RUNNER)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules["run_v9_candidate"] = MODULE
SPEC.loader.exec_module(MODULE)


def test_exact_clean_pair_uses_measured_residual():
    score, diag = MODULE.score_v9_clean_candidate(
        ("A1C", "B2D"),
        wt_target=1.0,
        clean_main={"A1C": 0.4, "B2D": -0.2},
        exact_pairs={("A1C", "B2D"): 0.7},
        relation_model={
            "bias": 100.0,
            "effect": {"A1C": 100.0, "B2D": 100.0},
        },
        calibration={"alpha": 100.0, "beta": 100.0},
        b2_fallback=-9.0,
    )

    assert np.isclose(score, 1.0 + 0.4 - 0.2 + 0.7)
    assert diag["exact_clean_pair_count"] == 1
    assert diag["calibrated_unseen_pair_count"] == 0


def test_unseen_pair_uses_affine_calibrated_relation():
    relation = {
        "bias": 0.5,
        "effect": {"A1C": 0.2, "B2D": -0.1},
    }
    calibration = {"alpha": 0.1, "beta": 0.25}

    score, diag = MODULE.score_v9_clean_candidate(
        ("A1C", "B2D"),
        wt_target=1.0,
        clean_main={"A1C": 0.4, "B2D": -0.2},
        exact_pairs={},
        relation_model=relation,
        calibration=calibration,
        b2_fallback=-9.0,
    )

    raw = 0.5 + 0.2 - 0.1
    pair = 0.1 + 0.25 * raw
    assert np.isclose(score, 1.0 + 0.4 - 0.2 + pair)
    assert diag["exact_clean_pair_count"] == 0
    assert diag["calibrated_unseen_pair_count"] == 1


def test_missing_singleton_falls_back_to_entire_b2_score():
    score, diag = MODULE.score_v9_clean_candidate(
        ("A1C", "X9Y"),
        wt_target=1.0,
        clean_main={"A1C": 0.4},
        exact_pairs={},
        relation_model={
            "bias": 0.0,
            "effect": {"A1C": 1.0},
        },
        calibration={"alpha": 0.0, "beta": 1.0},
        b2_fallback=3.25,
    )

    assert score == 3.25
    assert diag["used_b2_fallback"] is True
    assert diag["complete_clean_main"] is False


def test_second_order_expansion_sums_each_unordered_pair_once():
    relation = {
        "bias": 0.0,
        "effect": {
            "A1C": 1.0,
            "B2D": 2.0,
            "C3E": 3.0,
        },
    }
    calibration = {"alpha": 0.0, "beta": 1.0}

    score, diag = MODULE.score_v9_clean_candidate(
        ("A1C", "B2D", "C3E"),
        wt_target=0.0,
        clean_main={
            "A1C": 0.1,
            "B2D": 0.2,
            "C3E": 0.3,
        },
        exact_pairs={},
        relation_model=relation,
        calibration=calibration,
        b2_fallback=-9.0,
    )

    pair_sum = (1.0 + 2.0) + (1.0 + 3.0) + (2.0 + 3.0)
    assert np.isclose(score, 0.1 + 0.2 + 0.3 + pair_sum)
    assert diag["calibrated_unseen_pair_count"] == 3
    assert diag["unresolved_pair_count"] == 0
