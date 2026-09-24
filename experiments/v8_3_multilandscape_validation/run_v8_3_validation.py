import argparse
import importlib.util
import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).resolve()
V81_PATH = (
    HERE.parents[1]
    / "v8_1_crossfit_higher_order_dev"
    / "run_v8_1_crossfit_dev.py"
)

spec = importlib.util.spec_from_file_location("nabu_v81", V81_PATH)
if spec is None or spec.loader is None:
    raise RuntimeError(f"Could not import V8.1 runner from {V81_PATH}")
v81 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v81)


def percentile(values, x):
    values = np.asarray(values, dtype=float)
    return float(
        (np.sum(values < x) + 0.5 * np.sum(values == x))
        / len(values)
    )


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


def rank_preserving_elite_rerank(frame, fraction=0.20):
    n = len(frame)
    k = max(1, int(np.ceil(n * fraction)))
    b3_order = frame.sort_values(
        ["B3_RAW_PAIR", "candidate_id"],
        ascending=[False, True],
    ).index.tolist()
    elite = b3_order[:k]
    rest = b3_order[k:]
    elite_sorted = (
        frame.loc[elite]
        .sort_values(
            [
                "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER",
                "candidate_id",
            ],
            ascending=[False, True],
        )
        .index.tolist()
    )
    order = elite_sorted + rest
    score = pd.Series(index=frame.index, dtype=float)
    for rank, idx in enumerate(order):
        score.loc[idx] = float(n - rank)
    return score


def replacement(frame, arm, baseline, truth_column, k=50):
    arm_ids = set(
        frame.sort_values(
            [arm, "candidate_id"],
            ascending=[False, True],
        ).head(k)["candidate_id"]
    )
    base_ids = set(
        frame.sort_values(
            [baseline, "candidate_id"],
            ascending=[False, True],
        ).head(k)["candidate_id"]
    )
    added = arm_ids - base_ids
    removed = base_ids - arm_ids
    truth = frame.set_index("candidate_id")[truth_column].astype(float)
    added_mean = float(truth.loc[list(added)].mean()) if added else None
    removed_mean = float(truth.loc[list(removed)].mean()) if removed else None
    return {
        "overlap": int(len(arm_ids & base_ids)),
        "added_count": int(len(added)),
        "removed_count": int(len(removed)),
        "added_mean_true": added_mean,
        "removed_mean_true": removed_mean,
        "replacement_gain_true": (
            float(added_mean - removed_mean)
            if added_mean is not None and removed_mean is not None
            else None
        ),
    }


def run(csv_path, out_path):
    data = pd.read_csv(csv_path)
    mutation_column = v81.find_column(
        data.columns,
        v81.MUTATION_COLUMNS,
    )
    truth_column = v81.find_column(
        data.columns,
        v81.FITNESS_COLUMNS,
    )
    if mutation_column is None or truth_column is None:
        raise RuntimeError(
            f"Need mutation/fitness columns; found {list(data.columns)}"
        )

    data = data[[mutation_column, truth_column]].dropna().copy()
    data["mutation_set"] = data[mutation_column].map(v81.parse_mutations)
    data = data[
        data["mutation_set"].map(
            lambda x: 3 <= len(x) <= 5 and "*" not in "".join(x)
        )
    ].copy()
    data["candidate_id"] = data["mutation_set"].map(
        lambda x: ":".join(x)
    )
    data["bucket"] = data["candidate_id"].map(
        lambda x: v81.fnv1a32(x) % 10
    )

    visible = data[data["bucket"] <= 6].copy()
    hidden = data[data["bucket"] >= 7].copy()

    print(
        f"[V8.3] {Path(csv_path).name} "
        f"rows={len(data)} visible={len(visible)} hidden={len(hidden)}"
    )

    model = v81.fit_crossfitted_hierarchy(
        visible["mutation_set"].tolist(),
        visible[truth_column].to_numpy(dtype=float),
        visible["candidate_id"].tolist(),
    )

    eligible = hidden[
        hidden["mutation_set"].map(
            lambda mutation_set: all(
                mutation in model["base"]["main"]
                for mutation in mutation_set
            )
            and any(
                pair in model["base"]["pair"]
                for pair in combinations(mutation_set, 2)
            )
        )
    ].copy()

    scored = [
        v81.score_hierarchy(mutation_set, model)
        for mutation_set in eligible["mutation_set"]
    ]
    sf = pd.DataFrame(scored, index=eligible.index)
    for column in sf.columns:
        eligible[column] = sf[column]

    visible_labels = visible[truth_column].to_numpy(dtype=float)
    oof_b3 = float(
        spearmanr(model["oof_b3"], visible_labels).statistic
    )
    oof_b4 = float(
        spearmanr(model["oof_b4"], visible_labels).statistic
    )
    oof_delta = oof_b4 - oof_b3

    if oof_delta > 0:
        mode = "GLOBAL_HIGHER_ORDER"
        eligible["V8_3_ADAPTIVE_ROUTER"] = eligible[
            "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER"
        ].astype(float)
    else:
        mode = "RANK_PRESERVING_B3_TOP20_B5_RERANK"
        eligible["V8_3_ADAPTIVE_ROUTER"] = (
            rank_preserving_elite_rerank(eligible, 0.20)
        )

    arms = [
        "B3_RAW_PAIR",
        "B4_CROSSFIT_TRIPLET",
        "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER",
        "V8_3_ADAPTIVE_ROUTER",
    ]
    metrics = {
        arm: evaluate(eligible, arm, truth_column)
        for arm in arms
    }

    result = {
        "version": "NABU_V8_3_MULTILANDSCAPE_VALIDATION",
        "dataset": Path(csv_path).name,
        "rows": int(len(data)),
        "visible": int(len(visible)),
        "hidden": int(len(hidden)),
        "eligible_hidden": int(len(eligible)),
        "oof": {
            "B3_spearman": oof_b3,
            "B4_spearman": oof_b4,
            "B4_minus_B3": float(oof_delta),
            "router_mode": mode,
        },
        "metrics": metrics,
        "V8_3_vs_B3": {
            key: float(
                metrics["V8_3_ADAPTIVE_ROUTER"][key]
                - metrics["B3_RAW_PAIR"][key]
            )
            for key in [
                "spearman",
                "top1_percentile",
                "top5_mean_percentile",
                "top10_mean_percentile",
                "top50_mean_percentile",
                "top50_mean_true_score",
                "normalized_top1_regret",
            ]
        },
        "top50_replacement_vs_B3": replacement(
            eligible,
            "V8_3_ADAPTIVE_ROUTER",
            "B3_RAW_PAIR",
            truth_column,
            50,
        ),
        "model_diagnostics": v81.model_diagnostics(model),
    }

    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2), encoding="utf-8")

    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("csv")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    run(args.csv, args.out)
