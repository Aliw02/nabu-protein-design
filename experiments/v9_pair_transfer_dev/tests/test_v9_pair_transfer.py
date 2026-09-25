from __future__ import annotations

import importlib.util
from pathlib import Path
from types import SimpleNamespace

import numpy as np

MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "run_v9_pair_transfer_ired.py"
)
spec = importlib.util.spec_from_file_location("v9_pair_transfer", MODULE_PATH)
v9 = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(v9)


def test_factorization_recovers_missing_edge_from_node_propensities():
    mutations = [
        ("A1C", "B2D"),
        ("A1C", "C3E"),
        ("A1C", "D4F"),
        ("B2D", "C3E"),
        ("C3E", "D4F"),
    ]
    latent = {
        "A1C": 0.1,
        "B2D": 0.2,
        "C3E": 0.3,
        "D4F": 0.4,
    }
    bias = 1.0
    labels = np.array(
        [
            bias + latent[left] + latent[right]
            for left, right in mutations
        ],
        dtype=float,
    )
    fake = SimpleNamespace(
        model={"oof_b2": np.zeros(len(mutations), dtype=float)}
    )

    factor = v9.fit_pair_factorization(
        fake,
        mutations,
        labels,
        permute_targets=False,
    )
    predicted = v9.transferred_pair_effect(
        ("B2D", "D4F"),
        factor,
    )

    expected_raw = bias + latent["B2D"] + latent["D4F"]
    confidence = (
        factor["confidence"]["B2D"]
        * factor["confidence"]["D4F"]
    ) ** 0.5
    assert predicted is not None
    assert np.isclose(predicted, expected_raw * confidence)


def test_v9_preserves_exact_pair_and_transfers_only_missing_pair():
    fake = SimpleNamespace(
        model={
            "base": {
                "global_mean": 0.0,
                "main": {
                    "A1C": {"effect": 1.0},
                    "B2D": {"effect": 1.0},
                    "C3E": {"effect": 1.0},
                },
                "pair": {
                    ("A1C", "B2D"): {
                        "mean": 2.0,
                        "confidence": 0.5,
                    }
                },
            }
        }
    )
    factor = {
        "bias": 0.0,
        "propensity": {
            "A1C": 0.5,
            "B2D": 0.2,
            "C3E": 0.5,
        },
        "confidence": {
            "A1C": 1.0,
            "B2D": 1.0,
            "C3E": 1.0,
        },
    }

    exact_score, exact_diag = v9.score_v9(
        ("A1C", "B2D"),
        fake,
        factor,
    )
    transferred_score, transferred_diag = v9.score_v9(
        ("A1C", "C3E"),
        fake,
        factor,
    )

    assert np.isclose(exact_score, 3.0)
    assert exact_diag["exact_pair_count"] == 1
    assert exact_diag["transferred_pair_count"] == 0

    assert np.isclose(transferred_score, 3.0)
    assert transferred_diag["exact_pair_count"] == 0
    assert transferred_diag["transferred_pair_count"] == 1


def test_unseen_main_identity_keeps_neutral_abstention():
    fake = SimpleNamespace(
        model={
            "base": {
                "global_mean": 2.5,
                "main": {"A1C": {"effect": 1.0}},
                "pair": {},
            }
        }
    )
    factor = {
        "bias": 0.0,
        "propensity": {"A1C": 0.0},
        "confidence": {"A1C": 1.0},
    }

    score, diag = v9.score_v9(
        ("A1C", "B2D"),
        fake,
        factor,
    )

    assert score == 2.5
    assert diag["all_main_supported"] is False
