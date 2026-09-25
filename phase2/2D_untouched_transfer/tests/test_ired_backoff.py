from __future__ import annotations

from types import SimpleNamespace
from pathlib import Path
import sys

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from run_ired_blind import predict_with_frozen_backoff


class FakeModel:
    def __init__(self):
        self.model = {
            "base": {
                "global_mean": 10.0,
                "main": {
                    "A1C": {"effect": 2.0},
                    "A2C": {"effect": -1.0},
                    "A3C": {"effect": 3.0},
                },
                "pair": {
                    ("A1C", "A2C"): {
                        "mean": 1.0,
                        "confidence": 0.5,
                    }
                },
            }
        }

    def score_candidates(self, mutation_sets, candidate_ids):
        return pd.DataFrame(
            {
                "candidate_id": candidate_ids,
                "scoreable": [True] * len(candidate_ids),
                "V8_3_ADAPTIVE_ROUTER": [42.0] * len(candidate_ids),
            }
        )


def test_backoff_uses_deepest_supported_frozen_level():
    model = FakeModel()
    mutations = [
        ("A1C", "A2C"),
        ("A1C", "A3C"),
        ("A1C", "A4C"),
    ]
    result = predict_with_frozen_backoff(
        model,
        mutations,
        ["router", "b2", "abstain"],
    )

    assert result["prediction_source"].tolist() == [
        "V8_3_ADAPTIVE_ROUTER",
        "B2_ADDITIVE",
        "GLOBAL_MEAN_ABSTENTION",
    ]
    assert result["prediction"].tolist() == [
        42.0,
        15.0,
        10.0,
    ]


def test_backoff_is_finite_and_identity_only():
    model = FakeModel()
    result = predict_with_frozen_backoff(
        model,
        [
            ("A1C", "A2C"),
            ("A1C", "A3C"),
            ("A9C",),
        ],
        ["a", "b", "c"],
    )
    assert np.isfinite(
        result["prediction"].to_numpy(dtype=float)
    ).all()
    assert "target" not in result.columns
