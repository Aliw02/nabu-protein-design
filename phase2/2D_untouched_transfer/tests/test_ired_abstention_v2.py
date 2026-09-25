from __future__ import annotations

from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from run_ired_blind_v2 import apply_abstention_v2


def test_abstention_v2_preserves_scoreable_and_ties_unscoreable_below():
    router = np.array([2.5, np.nan, 1.25, np.nan], dtype=float)
    scoreable = np.array([True, False, True, False], dtype=bool)

    prediction, floor = apply_abstention_v2(router, scoreable)

    assert floor == 0.25
    assert prediction[0] == 2.5
    assert prediction[2] == 1.25
    assert prediction[1] == floor
    assert prediction[3] == floor
    assert floor < prediction[scoreable].min()


def test_abstention_v2_has_no_b2_or_target_input():
    router = np.array([5.0, np.nan, 3.0], dtype=float)
    scoreable = np.array([True, False, True], dtype=bool)

    prediction, _ = apply_abstention_v2(router, scoreable)

    assert np.isfinite(prediction).all()
    assert prediction.tolist() == [5.0, 2.0, 3.0]
