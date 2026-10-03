"""
Canonical evaluation helpers used by the frozen V8.3 validation path.
"""

from __future__ import annotations

import numpy as np
from scipy.stats import spearmanr


def percentile(values, x):
    values = np.asarray(values, dtype=float)
    return float(
        (np.sum(values < x) + 0.5 * np.sum(values == x))
        / len(values)
    )


def evaluate_v83(frame, score_column, truth_column):
    required = {"candidate_id", score_column, truth_column}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    if len(frame) == 0:
        raise ValueError("Cannot evaluate an empty candidate frame.")

    truth = frame[truth_column].to_numpy(dtype=float)
    scores = frame[score_column].to_numpy(dtype=float)
    if not np.isfinite(truth).all() or not np.isfinite(scores).all():
        raise ValueError("Truth and score columns must contain only finite values.")

    ranked = frame.sort_values(
        [score_column, "candidate_id"],
        ascending=[False, True],
    )
    top1_truth = float(ranked.iloc[0][truth_column])
    top5 = ranked.head(5)[truth_column].to_numpy(dtype=float)
    top10 = ranked.head(10)[truth_column].to_numpy(dtype=float)
    top50 = ranked.head(50)[truth_column].to_numpy(dtype=float)
    cutoff1 = float(np.quantile(truth, 0.99))
    cutoff01 = float(np.quantile(truth, 0.999))
    best = float(np.max(truth))
    worst = float(np.min(truth))

    return {
        "spearman": float(
            spearmanr(frame[score_column], truth).statistic
        ),
        "top1_percentile": percentile(truth, top1_truth),
        "top5_mean_percentile": float(
            np.mean([percentile(truth, x) for x in top5])
        ),
        "top10_mean_percentile": float(
            np.mean([percentile(truth, x) for x in top10])
        ),
        "top50_mean_percentile": float(
            np.mean([percentile(truth, x) for x in top50])
        ),
        "top10_top1pct_hits": int(np.sum(top10 >= cutoff1)),
        "top50_top1pct_hits": int(np.sum(top50 >= cutoff1)),
        "top10_top01pct_hits": int(np.sum(top10 >= cutoff01)),
        "top50_top01pct_hits": int(np.sum(top50 >= cutoff01)),
        "top50_mean_true_score": float(np.mean(top50)),
        "normalized_top1_regret": (
            float((best - top1_truth) / (best - worst))
            if best > worst
            else 0.0
        ),
    }
