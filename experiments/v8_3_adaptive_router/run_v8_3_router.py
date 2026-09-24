import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

DATASETS = [
    "PHOT_CHLRE_Chen_2023",
    "GB1_Wu2016_NABU_V8_1",
]


def percentile(values, x):
    values = np.asarray(values, dtype=float)
    return float((np.sum(values < x) + 0.5 * np.sum(values == x)) / len(values))


def evaluate(frame, score_column, truth_column):
    truth = frame[truth_column].to_numpy(dtype=float)
    ranked = frame.sort_values(
        [score_column, "candidate_id"],
        ascending=[False, True],
    )
    top1_truth = float(ranked.iloc[0][truth_column])
    top5 = ranked.head(5)[truth_column].to_numpy(dtype=float)
    top10 = ranked.head(10)[truth_column].to_numpy(dtype=float)
    top50 = ranked.head(50)[truth_column].to_numpy(dtype=float)
    top1_cutoff = float(np.quantile(truth, 0.99))
    best = float(np.max(truth))
    worst = float(np.min(truth))

    return {
        "spearman": float(spearmanr(frame[score_column], truth).statistic),
        "top1_percentile": percentile(truth, top1_truth),
        "top5_mean_percentile": float(np.mean([percentile(truth, x) for x in top5])),
        "top10_mean_percentile": float(np.mean([percentile(truth, x) for x in top10])),
        "top50_mean_percentile": float(np.mean([percentile(truth, x) for x in top50])),
        "top10_top1pct_hits": int(np.sum(top10 >= top1_cutoff)),
        "top50_top1pct_hits": int(np.sum(top50 >= top1_cutoff)),
        "top50_top1pct_enrichment": float(np.mean(top50 >= top1_cutoff) / 0.01),
        "top50_mean_true_score": float(np.mean(top50)),
        "normalized_top1_regret": (
            float((best - top1_truth) / (best - worst)) if best > worst else 0.0
        ),
    }


def infer_truth_column(frame):
    blocked = {
        "candidate_id",
        "bucket",
        "mutation_count",
        "B2_ADDITIVE",
        "B3_RAW_PAIR",
        "B4_CROSSFIT_TRIPLET",
        "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER",
        "triplet_delta",
        "quartet_delta",
        "triplet_supported",
        "quartet_supported",
        "triplet_confidence",
        "quartet_confidence",
        "GLOBAL_PERMUTED_B5",
        "COUNT_STRATIFIED_PERMUTED_B5",
    }
    blocked.update(c for c in frame.columns if c.startswith("T_"))

    preferred = ["DMS_score", "score", "fitness", "mean", "Fitness"]
    for column in preferred:
        if column in frame.columns:
            return column

    numeric = []
    for column in frame.columns:
        if column in blocked:
            continue
        if pd.api.types.is_numeric_dtype(frame[column]):
            numeric.append(column)

    if len(numeric) != 1:
        raise RuntimeError(
            f"Could not infer truth column. Numeric candidates: {numeric}"
        )
    return numeric[0]


def rank_preserving_elite_rerank(frame, fraction=0.20):
    n = len(frame)
    k = max(1, int(np.ceil(n * fraction)))

    b3_order = frame.sort_values(
        ["B3_RAW_PAIR", "candidate_id"],
        ascending=[False, True],
    ).index.tolist()

    elite_idx = b3_order[:k]
    rest_idx = b3_order[k:]

    elite_sorted = (
        frame.loc[elite_idx]
        .sort_values(
            ["B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER", "candidate_id"],
            ascending=[False, True],
        )
        .index.tolist()
    )

    final_order = elite_sorted + rest_idx

    # Score is rank-only by design. It guarantees:
    # 1) no candidate outside the B3 elite region crosses into it,
    # 2) B5 can only reorder candidates inside the frozen B3 elite region,
    # 3) outside-region B3 ordering remains exactly unchanged.
    score = pd.Series(index=frame.index, dtype=float)
    for rank, idx in enumerate(final_order):
        score.loc[idx] = float(n - rank)

    return score, {
        "elite_fraction": float(fraction),
        "elite_count": int(k),
        "boundary_preserved": True,
        "outside_B3_order_preserved": True,
    }


def replacement_diagnostic(frame, arm, baseline, truth_column, k=50):
    arm_ids = set(
        frame.sort_values([arm, "candidate_id"], ascending=[False, True])
        .head(k)["candidate_id"]
    )
    base_ids = set(
        frame.sort_values([baseline, "candidate_id"], ascending=[False, True])
        .head(k)["candidate_id"]
    )
    added = arm_ids - base_ids
    displaced = base_ids - arm_ids
    truth = frame.set_index("candidate_id")[truth_column].astype(float)

    added_mean = float(truth.loc[list(added)].mean()) if added else None
    displaced_mean = float(truth.loc[list(displaced)].mean()) if displaced else None

    return {
        "overlap": int(len(arm_ids & base_ids)),
        "added_count": int(len(added)),
        "displaced_count": int(len(displaced)),
        "added_mean_true": added_mean,
        "displaced_mean_true": displaced_mean,
        "replacement_gain_true": (
            float(added_mean - displaced_mean)
            if added_mean is not None and displaced_mean is not None
            else None
        ),
    }


def run_dataset(root, dataset):
    ddir = root / dataset
    scores_path = ddir / "V8_2_CANDIDATE_SCORES.csv"
    results_path = ddir / "V8_2_DATASET_RESULTS.json"

    frame = pd.read_csv(scores_path)
    prior = json.loads(results_path.read_text(encoding="utf-8"))
    truth_column = infer_truth_column(frame)

    oof_delta = float(
        prior["visible_oof_landscape_diagnostic"]["B4_minus_B3"]
    )

    if oof_delta > 0:
        frame["V8_3_ADAPTIVE_ROUTER"] = frame[
            "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER"
        ].astype(float)
        mode = "GLOBAL_HIGHER_ORDER"
        router_details = {
            "rule": "OOF_B4_MINUS_B3_POSITIVE",
            "boundary_preserved": None,
        }
    else:
        score, router_details = rank_preserving_elite_rerank(
            frame,
            fraction=0.20,
        )
        frame["V8_3_ADAPTIVE_ROUTER"] = score
        mode = "RANK_PRESERVING_B3_TOP20_B5_RERANK"

    arms = [
        "B3_RAW_PAIR",
        "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER",
        "T_OOF_OBJECTIVE_GATE",
        "V8_3_ADAPTIVE_ROUTER",
    ]
    metrics = {
        arm: evaluate(frame, arm, truth_column)
        for arm in arms
    }

    b3 = metrics["B3_RAW_PAIR"]
    v82 = metrics["T_OOF_OBJECTIVE_GATE"]
    v83 = metrics["V8_3_ADAPTIVE_ROUTER"]

    result = {
        "dataset": dataset,
        "truth_column": truth_column,
        "visible_oof_B4_minus_B3": oof_delta,
        "router_mode": mode,
        "router_details": router_details,
        "metrics": metrics,
        "V8_3_vs_B3": {
            key: float(v83[key] - b3[key])
            for key in [
                "spearman",
                "top1_percentile",
                "top5_mean_percentile",
                "top10_mean_percentile",
                "top50_mean_percentile",
                "top50_top1pct_enrichment",
                "top50_mean_true_score",
                "normalized_top1_regret",
            ]
        },
        "V8_3_vs_V8_2_GATE": {
            key: float(v83[key] - v82[key])
            for key in [
                "spearman",
                "top1_percentile",
                "top5_mean_percentile",
                "top10_mean_percentile",
                "top50_mean_percentile",
                "top50_top1pct_enrichment",
                "top50_mean_true_score",
                "normalized_top1_regret",
            ]
        },
        "top50_replacement_vs_B3": replacement_diagnostic(
            frame,
            "V8_3_ADAPTIVE_ROUTER",
            "B3_RAW_PAIR",
            truth_column,
            50,
        ),
        "top50_replacement_vs_V8_2_GATE": replacement_diagnostic(
            frame,
            "V8_3_ADAPTIVE_ROUTER",
            "T_OOF_OBJECTIVE_GATE",
            truth_column,
            50,
        ),
    }

    return result


def main(root_path, out_path):
    root = Path(root_path)
    out = Path(out_path)
    out.mkdir(parents=True, exist_ok=True)

    results = [run_dataset(root, dataset) for dataset in DATASETS]

    summary_rows = []
    for item in results:
        v83 = item["metrics"]["V8_3_ADAPTIVE_ROUTER"]
        b3 = item["metrics"]["B3_RAW_PAIR"]
        v82 = item["metrics"]["T_OOF_OBJECTIVE_GATE"]
        summary_rows.append(
            {
                "dataset": item["dataset"],
                "router_mode": item["router_mode"],
                "oof_delta": item["visible_oof_B4_minus_B3"],
                "B3_spearman": b3["spearman"],
                "V8_2_gate_spearman": v82["spearman"],
                "V8_3_spearman": v83["spearman"],
                "V8_3_delta_spearman_vs_B3": v83["spearman"] - b3["spearman"],
                "V8_3_delta_spearman_vs_V8_2": v83["spearman"] - v82["spearman"],
                "B3_top50_hits": b3["top50_top1pct_hits"],
                "V8_2_top50_hits": v82["top50_top1pct_hits"],
                "V8_3_top50_hits": v83["top50_top1pct_hits"],
                "B3_top50_true": b3["top50_mean_true_score"],
                "V8_2_top50_true": v82["top50_mean_true_score"],
                "V8_3_top50_true": v83["top50_mean_true_score"],
            }
        )

    summary = pd.DataFrame(summary_rows)

    final = {
        "version": "NABU_V8_3_ADAPTIVE_ROUTER",
        "status": "DEVELOPMENT_CANDIDATE",
        "architecture": (
            "B3 pairwise backbone + cross-fitted higher-order. "
            "Visible OOF selects global B5 when higher-order improves global OOF; "
            "otherwise B3 freezes top-20% membership and B5 reranks only inside "
            "that region with strict boundary preservation."
        ),
        "datasets": results,
        "aggregate": {
            "global_rank_gain_vs_B3_datasets": int(
                np.sum(summary["V8_3_delta_spearman_vs_B3"] > 0)
            ),
            "global_rank_gain_vs_V8_2_datasets": int(
                np.sum(summary["V8_3_delta_spearman_vs_V8_2"] > 0)
            ),
            "mean_delta_spearman_vs_B3": float(
                summary["V8_3_delta_spearman_vs_B3"].mean()
            ),
            "mean_delta_spearman_vs_V8_2": float(
                summary["V8_3_delta_spearman_vs_V8_2"].mean()
            ),
            "top50_hits_not_worse_than_V8_2_datasets": int(
                np.sum(summary["V8_3_top50_hits"] >= summary["V8_2_top50_hits"])
            ),
            "top50_true_not_worse_than_V8_2_datasets": int(
                np.sum(summary["V8_3_top50_true"] >= summary["V8_2_top50_true"])
            ),
        },
        "promotion_gate": {
            "promote_if": [
                "V8_3 Spearman > B3 on both landscapes",
                "V8_3 Top-50 hits >= V8.2 gate on both landscapes",
                "V8_3 Top-50 mean true score >= V8.2 gate on both landscapes",
                "GB1 boundary-preserving router improves or preserves V8.2 gate Spearman",
            ]
        },
    }

    (out / "V8_3_RESULTS.json").write_text(
        json.dumps(final, indent=2),
        encoding="utf-8",
    )
    summary.to_csv(out / "V8_3_SUMMARY.csv", index=False)
    print(json.dumps(final, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root",
        default="nabu_v8_2_tournament_results",
    )
    parser.add_argument(
        "--out",
        default="nabu_v8_3_router_results",
    )
    args = parser.parse_args()
    main(args.root, args.out)
