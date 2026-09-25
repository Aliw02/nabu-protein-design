"""
Canonical NABU V8.3 dual-objective router.

The rules are packaged from the frozen Phase-1 implementation. No thresholds
or routing semantics are changed here.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def percentile(values, x):
    values = np.asarray(values, dtype=float)
    return float(
        (np.sum(values < x) + 0.5 * np.sum(values == x))
        / len(values)
    )


def oof_elite_diagnostic(candidate_ids, labels, model, fraction=0.20):
    work = pd.DataFrame(
        {
            "candidate_id": list(candidate_ids),
            "truth": np.asarray(labels, dtype=float),
            "B3_OOF": np.asarray(model["oof_b3"], dtype=float),
            "B4_OOF": np.asarray(model["oof_b4"], dtype=float),
        }
    )

    n = len(work)
    k_region = max(1, int(math.ceil(n * fraction)))
    region = (
        work.sort_values(
            ["B3_OOF", "candidate_id"],
            ascending=[False, True],
        )
        .head(k_region)
        .copy()
    )

    top_k = min(50, len(region))
    b3_top = region.sort_values(
        ["B3_OOF", "candidate_id"],
        ascending=[False, True],
    ).head(top_k)
    b4_top = region.sort_values(
        ["B4_OOF", "candidate_id"],
        ascending=[False, True],
    ).head(top_k)

    all_truth = work["truth"].to_numpy(dtype=float)
    cutoff1 = float(np.quantile(all_truth, 0.99))

    def summarize(rows):
        vals = rows["truth"].to_numpy(dtype=float)
        return {
            "top50_mean_percentile": float(
                np.mean([percentile(all_truth, x) for x in vals])
            ),
            "top50_top1pct_hits": int(np.sum(vals >= cutoff1)),
            "top50_mean_visible_score": float(np.mean(vals)),
        }

    b3 = summarize(b3_top)
    b4 = summarize(b4_top)
    allow = (
        b4["top50_mean_percentile"] > b3["top50_mean_percentile"]
        and b4["top50_top1pct_hits"] >= b3["top50_top1pct_hits"]
    )

    return {
        "region_fraction": float(fraction),
        "region_count": int(k_region),
        "B3": b3,
        "B4": b4,
        "delta_top50_mean_percentile": float(
            b4["top50_mean_percentile"]
            - b3["top50_mean_percentile"]
        ),
        "delta_top50_top1pct_hits": int(
            b4["top50_top1pct_hits"]
            - b3["top50_top1pct_hits"]
        ),
        "allow_elite_rerank": bool(allow),
    }


def rank_preserving_elite_rerank(
    frame,
    b3_col,
    b5_col,
    output_col,
    fraction=0.20,
):
    n = len(frame)
    k = max(1, int(math.ceil(n * fraction)))

    b3_order = frame.sort_values(
        [b3_col, "candidate_id"],
        ascending=[False, True],
    ).index.tolist()
    elite = b3_order[:k]
    rest = b3_order[k:]

    elite_sorted = (
        frame.loc[elite]
        .sort_values(
            [b5_col, "candidate_id"],
            ascending=[False, True],
        )
        .index.tolist()
    )

    order = elite_sorted + rest
    score = pd.Series(index=frame.index, dtype=float)
    for rank, idx in enumerate(order):
        score.loc[idx] = float(n - rank)
    frame[output_col] = score

    return {
        "elite_fraction": float(fraction),
        "elite_count": int(k),
        "boundary_preserved": True,
        "outside_B3_order_preserved": True,
    }


def apply_router(
    frame,
    model,
    visible_ids,
    visible_labels,
    b3_col,
    b5_col,
    output_col,
):
    b3_oof = float(
        spearmanr(model["oof_b3"], visible_labels).statistic
    )
    b4_oof = float(
        spearmanr(model["oof_b4"], visible_labels).statistic
    )
    global_delta = b4_oof - b3_oof

    elite = oof_elite_diagnostic(
        visible_ids,
        visible_labels,
        model,
        fraction=0.20,
    )

    if global_delta > 0:
        frame[output_col] = frame[b5_col].astype(float)
        mode = "GLOBAL_HIGHER_ORDER"
        details = {"boundary_preserved": None}
    elif elite["allow_elite_rerank"]:
        details = rank_preserving_elite_rerank(
            frame,
            b3_col,
            b5_col,
            output_col,
            fraction=0.20,
        )
        mode = "RANK_PRESERVING_B3_TOP20_B5_RERANK"
    else:
        frame[output_col] = frame[b3_col].astype(float)
        mode = "B3_PROTECTED_NO_HIGHER_ORDER"
        details = {
            "boundary_preserved": True,
            "reason": (
                "Visible OOF did not support global or elite "
                "higher-order use."
            ),
        }

    return {
        "B3_oof_spearman": b3_oof,
        "B4_oof_spearman": b4_oof,
        "B4_minus_B3": float(global_delta),
        "elite_diagnostic": elite,
        "mode": mode,
        "details": details,
    }
