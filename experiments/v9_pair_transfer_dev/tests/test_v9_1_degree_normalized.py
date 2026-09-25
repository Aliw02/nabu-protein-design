from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[3]


def load_module(name: str, relative_path: str):
    path = ROOT / relative_path
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


V90 = load_module(
    "run_v9_pair_transfer_ired",
    "experiments/v9_pair_transfer_dev/run_v9_pair_transfer_ired.py",
)
V91 = load_module(
    "run_v9_1_degree_normalized",
    "experiments/v9_pair_transfer_dev/run_v9_1_degree_normalized.py",
)


class FakeModel:
    def __init__(self, base):
        self.model = {"base": base}


def test_single_transferred_edge_matches_v9_0():
    base = {
        "global_mean": 0.0,
        "main": {
            "A1C": {"effect": 0.0},
            "B2D": {"effect": 0.0},
        },
        "pair": {},
    }
    factor = {
        "bias": 0.5,
        "propensity": {
            "A1C": 1.0,
            "B2D": -0.25,
        },
        "confidence": {
            "A1C": 1.0,
            "B2D": 1.0,
        },
    }
    model = FakeModel(base)

    v90, _ = V90.score_v9(("A1C", "B2D"), model, factor)
    v91, diag = V91.score_v91(("A1C", "B2D"), model, factor)

    assert np.isclose(v91, v90)
    assert diag["transferred_pair_count"] == 1
    assert np.isclose(diag["transferred_mean_degree"], 1.0)


def test_complete_three_node_transfer_removes_degree_duplication():
    base = {
        "global_mean": 0.0,
        "main": {
            "A1C": {"effect": 0.0},
            "B2D": {"effect": 0.0},
            "C3E": {"effect": 0.0},
        },
        "pair": {},
    }
    factor = {
        "bias": 0.0,
        "propensity": {
            "A1C": 1.0,
            "B2D": 2.0,
            "C3E": 3.0,
        },
        "confidence": {
            "A1C": 1.0,
            "B2D": 1.0,
            "C3E": 1.0,
        },
    }
    model = FakeModel(base)

    v90, _ = V90.score_v9(("A1C", "B2D", "C3E"), model, factor)
    v91, diag = V91.score_v91(("A1C", "B2D", "C3E"), model, factor)

    assert np.isclose(v90, 12.0)
    assert np.isclose(diag["transferred_mean_degree"], 2.0)
    assert np.isclose(diag["raw_transferred_sum"], 12.0)
    assert np.isclose(diag["normalized_transferred_sum"], 6.0)
    assert np.isclose(v91, 6.0)


def test_exact_v83_pair_contribution_is_not_degree_normalized():
    base = {
        "global_mean": 0.0,
        "main": {
            "A1C": {"effect": 0.0},
            "B2D": {"effect": 0.0},
            "C3E": {"effect": 0.0},
        },
        "pair": {
            ("A1C", "B2D"): {
                "mean": 4.0,
                "confidence": 0.5,
            }
        },
    }
    factor = {
        "bias": 0.0,
        "propensity": {
            "A1C": 1.0,
            "B2D": 2.0,
            "C3E": 3.0,
        },
        "confidence": {
            "A1C": 1.0,
            "B2D": 1.0,
            "C3E": 1.0,
        },
    }
    model = FakeModel(base)

    score, diag = V91.score_v91(("A1C", "B2D", "C3E"), model, factor)

    # Exact A-B contributes +2 unchanged.
    # Transferred A-C and B-C contribute 4 + 5 = 9.
    # Their active graph has E=2, N=3, mean degree=4/3.
    expected_transfer = 9.0 / (4.0 / 3.0)
    assert np.isclose(score, 2.0 + expected_transfer)
    assert diag["exact_pair_count"] == 1
    assert diag["transferred_pair_count"] == 2


def test_unseen_main_still_abstains_to_global_mean():
    base = {
        "global_mean": 2.5,
        "main": {"A1C": {"effect": 0.0}},
        "pair": {},
    }
    factor = {
        "bias": 1.0,
        "propensity": {"A1C": 1.0},
        "confidence": {"A1C": 1.0},
    }
    model = FakeModel(base)

    score, diag = V91.score_v91(("A1C", "X9Y"), model, factor)

    assert score == 2.5
    assert diag["all_main_supported"] is False
    assert diag["transferred_pair_count"] == 0
