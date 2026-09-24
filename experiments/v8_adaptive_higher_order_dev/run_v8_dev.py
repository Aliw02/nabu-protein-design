import argparse
import json
import math
import re
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

MUTATION_COLUMNS = ["mutant", "mutation", "mutations", "variant"]
FITNESS_COLUMNS = ["DMS_score", "score", "fitness", "mean", "Fitness"]
SEED = 161


def fnv1a32(text: str) -> int:
    h = 2166136261
    for byte in text.encode("utf-8"):
        h = ((h ^ byte) * 16777619) & 0xFFFFFFFF
    return h


def parse_mutations(value):
    tokens = re.findall(r"[A-Z][0-9]+[A-Z*]", str(value))
    return tuple(sorted(tokens, key=lambda token: (int(token[1:-1]), token)))


def percentile(values, x):
    values = np.asarray(values, dtype=float)
    return float((np.sum(values < x) + 0.5 * np.sum(values == x)) / len(values))


def find_column(columns, candidates):
    return next((name for name in candidates if name in columns), None)


def build_subset_memory(mutation_sets, residuals, order, shrink_lambda, min_support):
    keys = []
    vals = []
    for mutation_set, residual in zip(mutation_sets, residuals):
        if len(mutation_set) < order:
            continue
        for subset in combinations(mutation_set, order):
            keys.append(subset)
            vals.append(float(residual))
    if not keys:
        return {}
    frame = pd.DataFrame({"key": pd.Series(keys, dtype="object"), "value": vals})
    grouped = frame.groupby("key", sort=False)["value"].agg(["mean", "count"])
    memory = {}
    for key, row in grouped.iterrows():
        support = int(row["count"])
        if support < min_support:
            continue
        confidence = support / (support + shrink_lambda)
        memory[key] = {
            "mean": float(row["mean"]),
            "support": support,
            "confidence": float(confidence),
        }
    return memory


def build_main_memory(mutation_sets, labels, global_mean, shrink_lambda):
    keys = []
    vals = []
    for mutation_set, label in zip(mutation_sets, labels):
        for mutation in mutation_set:
            keys.append(mutation)
            vals.append(float(label))
    frame = pd.DataFrame({"key": keys, "value": vals})
    grouped = frame.groupby("key", sort=False)["value"].agg(["mean", "count"])
    memory = {}
    for key, row in grouped.iterrows():
        support = int(row["count"])
        confidence = support / (support + shrink_lambda)
        effect = (float(row["mean"]) - global_mean) * confidence
        memory[key] = {
            "effect": float(effect),
            "support": support,
            "confidence": float(confidence),
        }
    return memory


def additive_score(mutation_set, global_mean, main_memory):
    return global_mean + sum(
        main_memory.get(mutation, {}).get("effect", 0.0) for mutation in mutation_set
    )


def pair_score(mutation_set, base_score, pair_memory):
    pair_effect = 0.0
    for pair in combinations(mutation_set, 2):
        entry = pair_memory.get(pair)
        if entry is not None:
            pair_effect += entry["mean"] * entry["confidence"]
    return base_score + pair_effect


def adaptive_order_contribution(mutation_set, memory, order):
    possible = math.comb(len(mutation_set), order) if len(mutation_set) >= order else 0
    if possible == 0:
        return 0.0, 0, 0.0
    entries = [
        memory[subset]
        for subset in combinations(mutation_set, order)
        if subset in memory
    ]
    if not entries:
        return 0.0, 0, 0.0
    confidence_sum = sum(entry["confidence"] for entry in entries)
    if confidence_sum <= 0:
        return 0.0, len(entries), 0.0
    weighted_residual = (
        sum(entry["mean"] * entry["confidence"] for entry in entries)
        / confidence_sum
    )
    adaptive_confidence = min(1.0, confidence_sum / possible)
    return (
        float(weighted_residual * adaptive_confidence),
        len(entries),
        float(adaptive_confidence),
    )


def fit_hierarchy(mutation_sets, labels):
    labels = np.asarray(labels, dtype=float)
    global_mean = float(np.mean(labels))
    main_memory = build_main_memory(
        mutation_sets, labels, global_mean, shrink_lambda=1.0
    )

    b2_visible = np.array(
        [
            additive_score(mutation_set, global_mean, main_memory)
            for mutation_set in mutation_sets
        ],
        dtype=float,
    )
    pair_residual = labels - b2_visible
    pair_memory = build_subset_memory(
        mutation_sets,
        pair_residual,
        order=2,
        shrink_lambda=1.0,
        min_support=1,
    )
    b3_visible = np.array(
        [
            pair_score(mutation_set, base, pair_memory)
            for mutation_set, base in zip(mutation_sets, b2_visible)
        ],
        dtype=float,
    )

    triplet_residual = labels - b3_visible
    triplet_memory = build_subset_memory(
        mutation_sets,
        triplet_residual,
        order=3,
        shrink_lambda=2.0,
        min_support=2,
    )
    triplet_contrib = np.array(
        [
            adaptive_order_contribution(mutation_set, triplet_memory, 3)[0]
            for mutation_set in mutation_sets
        ],
        dtype=float,
    )
    b4_visible = b3_visible + triplet_contrib

    quartet_residual = labels - b4_visible
    quartet_memory = build_subset_memory(
        mutation_sets,
        quartet_residual,
        order=4,
        shrink_lambda=4.0,
        min_support=2,
    )

    return {
        "global_mean": global_mean,
        "main": main_memory,
        "pair": pair_memory,
        "triplet": triplet_memory,
        "quartet": quartet_memory,
    }


def score_hierarchy(mutation_set, model):
    b2 = additive_score(mutation_set, model["global_mean"], model["main"])
    b3 = pair_score(mutation_set, b2, model["pair"])

    triplet, triplet_supported, triplet_confidence = adaptive_order_contribution(
        mutation_set, model["triplet"], 3
    )
    b4 = b3 + triplet

    quartet, quartet_supported, quartet_confidence = adaptive_order_contribution(
        mutation_set, model["quartet"], 4
    )
    b5 = b4 + quartet

    return {
        "B2_ADDITIVE": float(b2),
        "B3_RAW_PAIR": float(b3),
        "B4_TRIPLET_RESIDUAL": float(b4),
        "B5_ADAPTIVE_HIGHER_ORDER": float(b5),
        "triplet_supported": int(triplet_supported),
        "triplet_confidence": float(triplet_confidence),
        "quartet_supported": int(quartet_supported),
        "quartet_confidence": float(quartet_confidence),
    }


def evaluate(frame, score_column, truth_column):
    truth = frame[truth_column].to_numpy(dtype=float)
    ranked = frame.sort_values(
        [score_column, "candidate_id"], ascending=[False, True]
    )
    top1_truth = float(ranked.iloc[0][truth_column])
    top5_truth = ranked.head(5)[truth_column].to_numpy(dtype=float)
    top50_truth = ranked.head(50)[truth_column].to_numpy(dtype=float)
    top1_cutoff = float(np.quantile(truth, 0.99))
    best = float(np.max(truth))
    worst = float(np.min(truth))
    correlation = spearmanr(
        frame[score_column], frame[truth_column]
    ).statistic
    return {
        "spearman": float(correlation),
        "top1_percentile": percentile(truth, top1_truth),
        "top5_mean_percentile": float(
            np.mean([percentile(truth, value) for value in top5_truth])
        ),
        "top50_top1pct_enrichment": float(
            np.mean(top50_truth >= top1_cutoff) / 0.01
        ),
        "normalized_top1_regret": (
            float((best - top1_truth) / (best - worst))
            if best > worst
            else 0.0
        ),
    }


def overlap_count(frame, left, right, k):
    left_ids = set(
        frame.sort_values(
            [left, "candidate_id"], ascending=[False, True]
        ).head(k)["candidate_id"]
    )
    right_ids = set(
        frame.sort_values(
            [right, "candidate_id"], ascending=[False, True]
        ).head(k)["candidate_id"]
    )
    return len(left_ids & right_ids)


def permute_within_mutation_count(labels, mutation_sets, seed):
    labels = np.asarray(labels, dtype=float)
    counts = np.array(
        [len(mutation_set) for mutation_set in mutation_sets], dtype=int
    )
    result = labels.copy()
    rng = np.random.default_rng(seed)
    for count in sorted(set(counts.tolist())):
        idx = np.where(counts == count)[0]
        result[idx] = rng.permutation(labels[idx])
    return result


def memory_diagnostics(model):
    def summarize(memory):
        supports = np.array(
            [entry["support"] for entry in memory.values()], dtype=float
        )
        if len(supports) == 0:
            return {
                "entries": 0,
                "support_min": 0,
                "support_median": 0.0,
                "support_max": 0,
            }
        return {
            "entries": int(len(memory)),
            "support_min": int(np.min(supports)),
            "support_median": float(np.median(supports)),
            "support_max": int(np.max(supports)),
        }

    return {
        "main": summarize(model["main"]),
        "pair": summarize(model["pair"]),
        "triplet": summarize(model["triplet"]),
        "quartet": summarize(model["quartet"]),
    }


def main(csv_path, out_dir):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    data = pd.read_csv(csv_path)
    mutation_column = find_column(data.columns, MUTATION_COLUMNS)
    truth_column = find_column(data.columns, FITNESS_COLUMNS)
    if mutation_column is None or truth_column is None:
        raise SystemExit(
            f"Need mutation and fitness columns; found {list(data.columns)}"
        )

    data = data[[mutation_column, truth_column]].dropna().copy()
    data["mutation_set"] = data[mutation_column].map(parse_mutations)
    data = data[
        data["mutation_set"].map(
            lambda x: 3 <= len(x) <= 5 and "*" not in "".join(x)
        )
    ].copy()
    data["candidate_id"] = data["mutation_set"].map(lambda x: ":".join(x))
    data["bucket"] = data["candidate_id"].map(lambda x: fnv1a32(x) % 10)
    data["mutation_count"] = data["mutation_set"].map(len)

    visible = data[data["bucket"] <= 6].copy()
    hidden = data[data["bucket"] >= 7].copy()

    visible_sets = visible["mutation_set"].tolist()
    visible_labels = visible[truth_column].to_numpy(dtype=float)
    model = fit_hierarchy(visible_sets, visible_labels)

    eligible = hidden[
        hidden["mutation_set"].map(
            lambda mutation_set: all(
                mutation in model["main"] for mutation in mutation_set
            )
            and any(
                pair in model["pair"]
                for pair in combinations(mutation_set, 2)
            )
        )
    ].copy()

    scored = [
        score_hierarchy(mutation_set, model)
        for mutation_set in eligible["mutation_set"]
    ]
    scored_frame = pd.DataFrame(scored, index=eligible.index)
    for column in scored_frame.columns:
        eligible[column] = scored_frame[column]

    global_perm_labels = np.random.default_rng(SEED).permutation(
        visible_labels
    )
    global_perm_model = fit_hierarchy(visible_sets, global_perm_labels)
    eligible["GLOBAL_PERMUTED_B5"] = [
        score_hierarchy(
            mutation_set, global_perm_model
        )["B5_ADAPTIVE_HIGHER_ORDER"]
        for mutation_set in eligible["mutation_set"]
    ]

    count_perm_labels = permute_within_mutation_count(
        visible_labels, visible_sets, SEED
    )
    count_perm_model = fit_hierarchy(visible_sets, count_perm_labels)
    eligible["COUNT_STRATIFIED_PERMUTED_B5"] = [
        score_hierarchy(
            mutation_set, count_perm_model
        )["B5_ADAPTIVE_HIGHER_ORDER"]
        for mutation_set in eligible["mutation_set"]
    ]

    arms = [
        "B2_ADDITIVE",
        "B3_RAW_PAIR",
        "B4_TRIPLET_RESIDUAL",
        "B5_ADAPTIVE_HIGHER_ORDER",
        "GLOBAL_PERMUTED_B5",
        "COUNT_STRATIFIED_PERMUTED_B5",
    ]
    metrics = {
        arm: evaluate(eligible, arm, truth_column)
        for arm in arms
    }

    by_mutation_count = {}
    for mutation_count, group in eligible.groupby("mutation_count"):
        if len(group) < 10:
            continue
        by_mutation_count[str(int(mutation_count))] = {
            arm: evaluate(group, arm, truth_column)
            for arm in arms[:4]
        }

    result = {
        "version": "NABU_V8_ADAPTIVE_HIGHER_ORDER_DEV",
        "status": "DEVELOPMENT_DIAGNOSTIC_NOT_FROZEN",
        "dataset": Path(csv_path).name,
        "rows": int(len(data)),
        "visible": int(len(visible)),
        "hidden": int(len(hidden)),
        "eligible_hidden": int(len(eligible)),
        "metrics": metrics,
        "deltas_vs_raw_pair": {
            "B4_minus_B3_spearman": float(
                metrics["B4_TRIPLET_RESIDUAL"]["spearman"]
                - metrics["B3_RAW_PAIR"]["spearman"]
            ),
            "B5_minus_B3_spearman": float(
                metrics["B5_ADAPTIVE_HIGHER_ORDER"]["spearman"]
                - metrics["B3_RAW_PAIR"]["spearman"]
            ),
            "B5_minus_B3_top5_mean_percentile": float(
                metrics["B5_ADAPTIVE_HIGHER_ORDER"]["top5_mean_percentile"]
                - metrics["B3_RAW_PAIR"]["top5_mean_percentile"]
            ),
            "B5_minus_B3_top50_enrichment": float(
                metrics["B5_ADAPTIVE_HIGHER_ORDER"][
                    "top50_top1pct_enrichment"
                ]
                - metrics["B3_RAW_PAIR"]["top50_top1pct_enrichment"]
            ),
        },
        "candidate_overlap_B3_vs_B5": {
            "top1": overlap_count(
                eligible, "B3_RAW_PAIR", "B5_ADAPTIVE_HIGHER_ORDER", 1
            ),
            "top5": overlap_count(
                eligible, "B3_RAW_PAIR", "B5_ADAPTIVE_HIGHER_ORDER", 5
            ),
            "top10": overlap_count(
                eligible, "B3_RAW_PAIR", "B5_ADAPTIVE_HIGHER_ORDER", 10
            ),
            "top20": overlap_count(
                eligible, "B3_RAW_PAIR", "B5_ADAPTIVE_HIGHER_ORDER", 20
            ),
            "top50": overlap_count(
                eligible, "B3_RAW_PAIR", "B5_ADAPTIVE_HIGHER_ORDER", 50
            ),
        },
        "higher_order_hidden_coverage": {
            "triplet_any_fraction": float(
                np.mean(eligible["triplet_supported"] > 0)
            ),
            "quartet_any_fraction": float(
                np.mean(eligible["quartet_supported"] > 0)
            ),
            "triplet_mean_confidence": float(
                np.mean(eligible["triplet_confidence"])
            ),
            "quartet_mean_confidence": float(
                np.mean(eligible["quartet_confidence"])
            ),
        },
        "memory_diagnostics": memory_diagnostics(model),
        "by_mutation_count": by_mutation_count,
        "interpretation_boundary": (
            "Development diagnostic on a previously revealed assay. "
            "Do not treat this run as fresh blind confirmation or a "
            "canonical breakpoint."
        ),
    }

    (out / "V8_DEV_RESULTS.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )

    export = eligible.drop(columns=["mutation_set"]).copy()
    export.to_csv(
        out / "V8_DEV_REVEALED_PREDICTIONS.csv", index=False
    )

    top_rows = []
    for arm in arms[:4]:
        ranked = eligible.sort_values(
            [arm, "candidate_id"], ascending=[False, True]
        ).head(50)
        for rank, row in enumerate(ranked.itertuples(), start=1):
            top_rows.append(
                {
                    "arm": arm,
                    "rank": rank,
                    "candidate_id": row.candidate_id,
                    "predicted_score": float(getattr(row, arm)),
                    "true_score": float(getattr(row, truth_column)),
                }
            )
    pd.DataFrame(top_rows).to_csv(
        out / "V8_DEV_TOP50.csv", index=False
    )

    (out / "V8_DEV_MEMORY_DIAGNOSTICS.json").write_text(
        json.dumps(memory_diagnostics(model), indent=2),
        encoding="utf-8",
    )

    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("csv")
    parser.add_argument("--out", default="nabu_v8_dev_results")
    args = parser.parse_args()
    main(args.csv, args.out)
