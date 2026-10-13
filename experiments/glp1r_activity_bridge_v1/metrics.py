"""Prespecified metrics and paired resampling, shared by freeze and evaluation."""
from __future__ import annotations

import numpy as np
from scipy.stats import spearmanr


def metrics(y, predictions, identifiers):
    y, predictions = np.asarray(y), np.asarray(predictions)
    assert len(y) == len(predictions) == len(identifiers) and np.isfinite(predictions).all()
    count = min(5, len(y))
    predicted_best = sorted(range(len(y)), key=lambda i: (-predictions[i], str(identifiers[i])))[:count]
    true_best = sorted(range(len(y)), key=lambda i: (-y[i], str(identifiers[i])))[:count]
    rank = None if np.ptp(y) == 0 or np.ptp(predictions) == 0 else float(spearmanr(y, predictions).statistic)
    return {"n": len(y), "coverage": 1.0,
            "rmse": float(np.sqrt(np.mean((predictions - y) ** 2))),
            "mae": float(np.mean(np.abs(predictions - y))), "spearman": rank,
            "top5_mean_measured_pEC50": float(y[predicted_best].mean()),
            "top5_recall": len(set(predicted_best).intersection(true_best)) / count,
            "top5_ids": [str(identifiers[i]) for i in predicted_best]}


def paired_bootstrap(y, candidate, baseline, groups, repeats, seed):
    delta = (np.asarray(candidate) - y) ** 2 - (np.asarray(baseline) - y) ** 2
    units = [np.flatnonzero(groups == group) for group in np.unique(groups)]
    rng = np.random.default_rng(seed)
    estimates = np.empty(repeats)
    for sample in range(repeats):
        selected = np.concatenate([units[i] for i in rng.integers(0, len(units), size=len(units))])
        estimates[sample] = delta[selected].mean()
    low, high = np.quantile(estimates, [0.025, 0.975])
    return {"candidate_minus_baseline_mse": float(delta.mean()), "ci95_low": float(low),
            "ci95_high": float(high), "resampling_units": len(units), "resamples": repeats,
            "upper_bound_negative": bool(high < 0)}
