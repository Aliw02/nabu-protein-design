from __future__ import annotations

import numpy as np
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from run_ired_blind_v2 import evaluate_complete_test


def test_complete_test_metrics_perfect_ranking():
    target = np.arange(100, dtype=float)
    prediction = target.copy()
    metrics = evaluate_complete_test(target, prediction)
    assert np.isclose(metrics["spearman"], 1.0)
    assert abs(metrics["ndcg"] - 1.0) < 1e-12
    assert metrics["top1_percent_recall"] == 1.0
    assert metrics["normalized_regret_top1pct"] == 0.0


def test_complete_test_metrics_do_not_subset():
    target = np.array([0.0, 1.0, 2.0, 3.0])
    prediction = np.array([3.0, 2.0, 1.0, 0.0])
    metrics = evaluate_complete_test(target, prediction)
    assert np.isclose(metrics["spearman"], -1.0)
    assert metrics["top1_percent_k"] == 1
