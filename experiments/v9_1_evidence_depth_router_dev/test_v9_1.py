from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from experiments.v9_1_evidence_depth_router_dev.run_v9_1_ired import (
    b2_b3_predictions,
)


def test_no_exact_pair_reduces_b3_to_b2():
    model = SimpleNamespace(
        model={
            "base": {
                "global_mean": 0.0,
                "main": {
                    "A1C": {"effect": 1.0},
                    "B2D": {"effect": 2.0},
                },
                "pair": {},
            }
        }
    )
    result = b2_b3_predictions(model, [("A1C", "B2D")])
    assert result["all_main"].tolist() == [True]
    assert result["exact_pair"].tolist() == [False]
    assert result["b2"].tolist() == [3.0]
    assert result["b3"].tolist() == [3.0]


def test_exact_pair_adds_only_frozen_pair_effect():
    model = SimpleNamespace(
        model={
            "base": {
                "global_mean": 0.0,
                "main": {
                    "A1C": {"effect": 1.0},
                    "B2D": {"effect": 2.0},
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
    result = b2_b3_predictions(model, [("A1C", "B2D")])
    assert result["b2"].tolist() == [3.0]
    assert result["b3"].tolist() == [5.0]
    assert result["exact_pair"].tolist() == [True]
    assert result["exact_pair_count"].tolist() == [1]


def test_unseen_main_uses_global_mean_for_both_depths():
    model = SimpleNamespace(
        model={
            "base": {
                "global_mean": 7.5,
                "main": {"A1C": {"effect": 1.0}},
                "pair": {},
            }
        }
    )
    result = b2_b3_predictions(model, [("A1C", "B2D")])
    assert result["all_main"].tolist() == [False]
    assert result["b2"].tolist() == [7.5]
    assert result["b3"].tolist() == [7.5]
    assert np.isfinite(result["b2"]).all()
    assert np.isfinite(result["b3"]).all()
