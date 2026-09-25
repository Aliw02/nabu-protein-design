from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from experiments.v9_pair_transfer_dev.run_v9_pair_transfer import (
    fit_factorized_pair_transfer,
    score_v9,
    transferred_pair_value,
)


def test_factorized_pair_transfer_is_deterministic_and_finite():
    mutation_sets = [
        ("A1C", "B2D"),
        ("A1C", "C3E"),
        ("B2D", "D4F"),
        ("C3E", "D4F"),
    ]
    labels = np.array([1.0, 2.0, 3.0, 4.0])
    oof_b2 = np.zeros(4, dtype=float)

    left = fit_factorized_pair_transfer(
        mutation_sets,
        labels,
        oof_b2,
        permute=False,
    )
    right = fit_factorized_pair_transfer(
        mutation_sets,
        labels,
        oof_b2,
        permute=False,
    )

    assert left["bias"] == right["bias"]
    assert left["node_effect"] == right["node_effect"]
    assert left["node_confidence"] == right["node_confidence"]

    unseen = transferred_pair_value(("A1C", "D4F"), left)
    assert unseen is not None
    assert np.isfinite(unseen)


def test_permuted_control_changes_factor_solution():
    mutation_sets = [
        ("A1C", "B2D"),
        ("A1C", "C3E"),
        ("A1C", "D4F"),
        ("B2D", "C3E"),
        ("B2D", "D4F"),
    ]
    labels = np.array([0.1, 1.7, -0.2, 3.4, 0.8])
    oof_b2 = np.zeros(5, dtype=float)

    real = fit_factorized_pair_transfer(
        mutation_sets,
        labels,
        oof_b2,
        permute=False,
    )
    control = fit_factorized_pair_transfer(
        mutation_sets,
        labels,
        oof_b2,
        permute=True,
    )

    real_vector = np.array(
        [real["node_effect"][key] for key in sorted(real["node_effect"])]
    )
    control_vector = np.array(
        [control["node_effect"][key] for key in sorted(control["node_effect"])]
    )
    assert not np.array_equal(real_vector, control_vector)


def test_exact_pair_memory_has_priority_over_transfer():
    frozen = SimpleNamespace(
        model={
            "base": {
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
        }
    )
    factor = {
        "bias": 100.0,
        "node_effect": {"A1C": 0.0, "B2D": 0.0, "C3E": 0.0},
        "node_confidence": {"A1C": 1.0, "B2D": 1.0, "C3E": 1.0},
    }

    prediction, diagnostic = score_v9(
        [("A1C", "B2D"), ("A1C", "C3E")],
        frozen,
        factor,
    )

    assert prediction[0] == 2.0
    assert diagnostic[0]["exact_pair_count"] == 1
    assert diagnostic[0]["transferred_pair_count"] == 0

    assert prediction[1] == 100.0
    assert diagnostic[1]["exact_pair_count"] == 0
    assert diagnostic[1]["transferred_pair_count"] == 1


def test_unseen_main_identity_abstains_to_global_mean():
    frozen = SimpleNamespace(
        model={
            "base": {
                "global_mean": 7.5,
                "main": {"A1C": {"effect": 1.0}},
                "pair": {},
            }
        }
    )
    factor = {
        "bias": 1.0,
        "node_effect": {"A1C": 1.0},
        "node_confidence": {"A1C": 0.5},
    }

    prediction, diagnostic = score_v9(
        [("A1C", "B2D")],
        frozen,
        factor,
    )

    assert prediction.tolist() == [7.5]
    assert diagnostic[0]["all_main_supported"] is False
    assert diagnostic[0]["transferred_pair_count"] == 0
