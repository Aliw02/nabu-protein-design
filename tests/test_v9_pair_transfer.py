from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np


RUNNER_PATH = (
    Path(__file__).resolve().parents[1]
    / "experiments"
    / "v9_pair_transfer_dev"
    / "run_v9_pair_transfer.py"
)
SPEC = importlib.util.spec_from_file_location("run_v9_pair_transfer", RUNNER_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_factorized_pair_model_recovers_unseen_edge_in_additive_pair_system():
    mutations = [
        ("A1C", "B2D"),
        ("A1C", "C3E"),
        ("A1C", "D4F"),
        ("B2D", "C3E"),
        ("B2D", "D4F"),
    ]
    latent = {
        "A1C": 0.4,
        "B2D": -0.2,
        "C3E": 0.7,
        "D4F": 0.1,
    }
    bias = 0.3
    residual = np.asarray(
        [
            bias + latent[left] + latent[right]
            for left, right in mutations
        ],
        dtype=float,
    )
    labels = residual.copy()
    oof_b2 = np.zeros(len(labels), dtype=float)

    model = MODULE.build_pair_factor_model(
        mutations,
        labels,
        oof_b2,
        permute_targets=False,
    )
    predicted = MODULE.transferred_pair_effect(("C3E", "D4F"), model)

    raw = (
        model["bias"]
        + model["node_effect"]["C3E"]
        + model["node_effect"]["D4F"]
    )
    conf = (
        model["node_support"]["C3E"]
        / (model["node_support"]["C3E"] + 1.0)
        * model["node_support"]["D4F"]
        / (model["node_support"]["D4F"] + 1.0)
    ) ** 0.5

    assert predicted is not None
    assert np.isclose(predicted, raw * conf)


def test_v9_keeps_exact_v83_pair_memory_instead_of_transfer_override():
    base = {
        "global_mean": 0.0,
        "main": {
            "A1C": {"effect": 0.0},
            "B2D": {"effect": 0.0},
        },
        "pair": {
            ("A1C", "B2D"): {
                "mean": 2.0,
                "confidence": 0.5,
            }
        },
    }
    factor = {
        "bias": 100.0,
        "node_effect": {"A1C": 100.0, "B2D": 100.0},
        "node_support": {"A1C": 10, "B2D": 10},
    }

    score, diag = MODULE.score_v9_candidate(
        ("A1C", "B2D"),
        base,
        factor,
    )

    assert np.isclose(score, 1.0)
    assert diag["exact_pairs"] == 1
    assert diag["transferred_pairs"] == 0
    assert diag["unresolved_pairs"] == 0


def test_unseen_main_identity_uses_global_mean_without_partial_extrapolation():
    base = {
        "global_mean": 3.5,
        "main": {"A1C": {"effect": 1.0}},
        "pair": {},
    }
    factor = {
        "bias": 1.0,
        "node_effect": {"A1C": 1.0},
        "node_support": {"A1C": 2},
    }

    score, diag = MODULE.score_v9_candidate(
        ("A1C", "X9Y"),
        base,
        factor,
    )

    assert score == 3.5
    assert diag["all_main_supported"] is False
    assert diag["exact_pairs"] == 0
    assert diag["transferred_pairs"] == 0


def test_permuted_factor_control_is_deterministic():
    mutations = [
        ("A1C", "B2D"),
        ("A1C", "C3E"),
        ("B2D", "C3E"),
        ("A1C", "D4F"),
    ]
    labels = np.asarray([0.1, 0.4, -0.2, 1.3], dtype=float)
    oof_b2 = np.zeros(len(labels), dtype=float)

    left = MODULE.build_pair_factor_model(
        mutations,
        labels,
        oof_b2,
        permute_targets=True,
    )
    right = MODULE.build_pair_factor_model(
        mutations,
        labels,
        oof_b2,
        permute_targets=True,
    )

    assert left["bias"] == right["bias"]
    assert left["node_effect"] == right["node_effect"]
    assert left["node_support"] == right["node_support"]
