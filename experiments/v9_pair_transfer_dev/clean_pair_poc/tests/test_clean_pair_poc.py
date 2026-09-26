from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import numpy as np


RUNNER = (
    Path(__file__).resolve().parents[1]
    / "run_clean_pair_poc.py"
)
SPEC = importlib.util.spec_from_file_location("run_clean_pair_poc", RUNNER)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules["run_clean_pair_poc"] = MODULE
SPEC.loader.exec_module(MODULE)


def test_clean_pair_identity_formula():
    wt = 1.0
    single_a = 1.4
    single_b = 0.7
    double_ab = 1.5

    main_a = single_a - wt
    main_b = single_b - wt
    additive = wt + main_a + main_b
    residual = double_ab - additive

    assert np.isclose(additive, single_a + single_b - wt)
    assert np.isclose(residual, double_ab - single_a - single_b + wt)
    assert np.isclose(residual, 0.4)


def test_fold_assignment_is_deterministic_and_bounded():
    candidate = "ACDEFGHIK"
    left = MODULE.assign_fold(candidate)
    right = MODULE.assign_fold(candidate)

    assert left == right
    assert 0 <= left < MODULE.N_FOLDS


def test_node_additive_fit_predicts_training_relation():
    pairs = [
        ("A1C", "B2D"),
        ("A1C", "C3E"),
        ("B2D", "C3E"),
        ("A1C", "D4F"),
        ("B2D", "D4F"),
        ("C3E", "D4F"),
    ]
    latent = {
        "A1C": 0.2,
        "B2D": -0.1,
        "C3E": 0.4,
        "D4F": 0.3,
    }
    bias = -0.2
    residual = np.asarray(
        [
            bias + latent[left] + latent[right]
            for left, right in pairs
        ],
        dtype=float,
    )

    model = MODULE.fit_node_additive(
        pairs,
        residual,
        permute=False,
        fold=0,
    )

    prediction = np.asarray(
        [
            MODULE.predict_node_additive(pair, model)
            for pair in pairs
        ],
        dtype=float,
    )

    assert np.allclose(prediction, residual)


def test_missing_node_returns_none():
    model = {
        "bias": 0.0,
        "effect": {
            "A1C": 0.1,
            "B2D": 0.2,
        },
    }

    assert MODULE.predict_node_additive(("A1C", "B2D"), model) is not None
    assert MODULE.predict_node_additive(("A1C", "X9Y"), model) is None


def test_permuted_control_is_deterministic_per_fold():
    pairs = [
        ("A1C", "B2D"),
        ("A1C", "C3E"),
        ("B2D", "C3E"),
        ("A1C", "D4F"),
    ]
    residual = np.asarray([0.1, 0.4, -0.2, 1.3], dtype=float)

    left = MODULE.fit_node_additive(
        pairs,
        residual,
        permute=True,
        fold=2,
    )
    right = MODULE.fit_node_additive(
        pairs,
        residual,
        permute=True,
        fold=2,
    )

    assert left["bias"] == right["bias"]
    assert left["effect"] == right["effect"]
