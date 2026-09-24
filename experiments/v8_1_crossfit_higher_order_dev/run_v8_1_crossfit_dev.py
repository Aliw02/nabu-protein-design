import argparse
import json
import math
import re
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

SEED = 161
N_FOLDS = 5
MUTATION_COLUMNS = ["mutant", "mutation", "mutations", "variant"]
FITNESS_COLUMNS = ["DMS_score", "score", "fitness", "mean", "Fitness"]


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


def crossfit_fold(candidate_id: str) -> int:
    return fnv1a32(candidate_id + "|NABU_V8_1_CF|161") % N_FOLDS


def build_main_memory(mutation_sets, labels, global_mean, shrink_lambda=1.0):
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
        memory[key] = {
            "effect": float((float(row["mean"]) - global_mean) * confidence),
            "support": support,
            "confidence": float(confidence),
        }
    return memory


def build_subset_memory(
    mutation_sets,
    residuals,
    order,
    shrink_lambda,
    min_support,
):
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

    frame = pd.DataFrame(
        {"key": pd.Series(keys, dtype="object"), "value": vals}
    )
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


def additive_score(mutation_set, global_mean, main_memory):
    return global_mean + sum(
        main_memory.get(mutation, {}).get("effect", 0.0)
        for mutation in mutation_set
    )


def pair_score(mutation_set, base_score, pair_memory):
    pair_effect = 0.0
    for pair in combinations(mutation_set, 2):
        entry = pair_memory.get(pair)
        if entry is not None:
            pair_effect += entry["mean"] * entry["confidence"]
    return float(base_score + pair_effect)


def adaptive_order_contribution(mutation_set, memory, order):
    if len(mutation_set) < order:
        return 0.0, 0, 0.0

    possible = math.comb(len(mutation_set), order)
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
    coverage_confidence = min(1.0, confidence_sum / possible)
    return (
        float(weighted_residual * coverage_confidence),
        int(len(entries)),
        float(coverage_confidence),
    )


def fit_base(mutation_sets, labels):
    labels = np.asarray(labels, dtype=float)
    global_mean = float(np.mean(labels))
    main_memory = build_main_memory(
        mutation_sets, labels, global_mean, shrink_lambda=1.0
    )
    b2 = np.array(
        [
            additive_score(mutation_set, global_mean, main_memory)
            for mutation_set in mutation_sets
        ],
        dtype=float,
    )
    residual = labels - b2
    pair_memory = build_subset_memory(
        mutation_sets,
        residual,
        order=2,
        shrink_lambda=1.0,
        min_support=1,
    )
    return {
        "global_mean": global_mean,
        "main": main_memory,
        "pair": pair_memory,
    }


def score_base(mutation_set, base_model):
    b2 = additive_score(
        mutation_set,
        base_model["global_mean"],
        base_model["main"],
    )
    b3 = pair_score(mutation_set, b2, base_model["pair"])
    return float(b2), float(b3)


def fit_crossfitted_hierarchy(
    mutation_sets,
    labels,
    candidate_ids,
):
    labels = np.asarray(labels, dtype=float)
    mutation_sets = list(mutation_sets)
    candidate_ids = list(candidate_ids)
    n = len(labels)

    folds = np.array(
        [crossfit_fold(candidate_id) for candidate_id in candidate_ids],
        dtype=int,
    )

    if set(folds.tolist()) != set(range(N_FOLDS)):
        raise RuntimeError(
            "Cross-fit fold assignment did not populate all five folds."
        )

    oof_b2 = np.full(n, np.nan, dtype=float)
    oof_b3 = np.full(n, np.nan, dtype=float)

    # First cross-fit: unbiased B3 residual target for every visible row.
    for fold in range(N_FOLDS):
        train_idx = np.where(folds != fold)[0]
        hold_idx = np.where(folds == fold)[0]
        model = fit_base(
            [mutation_sets[i] for i in train_idx],
            labels[train_idx],
        )
        for i in hold_idx:
            b2, b3 = score_base(mutation_sets[i], model)
            oof_b2[i] = b2
            oof_b3[i] = b3

    if np.isnan(oof_b3).any():
        raise RuntimeError("Missing B3 OOF predictions.")

    residual3_oof = labels - oof_b3

    final_triplet_memory = build_subset_memory(
        mutation_sets,
        residual3_oof,
        order=3,
        shrink_lambda=2.0,
        min_support=2,
    )

    # Strict outer cross-fit B4 prediction.
    # For each outer fold, triplet memory is learned only from the other folds.
    # The residual targets used to build that triplet memory are themselves
    # generated by an inner cross-fit that excludes both the outer fold and
    # the row's own inner fold.
    oof_b4 = np.full(n, np.nan, dtype=float)

    for outer_fold in range(N_FOLDS):
        outer_train_mask = folds != outer_fold
        outer_hold_idx = np.where(folds == outer_fold)[0]
        outer_train_idx = np.where(outer_train_mask)[0]

        inner_residual3 = np.full(n, np.nan, dtype=float)

        for inner_fold in range(N_FOLDS):
            if inner_fold == outer_fold:
                continue

            inner_hold_idx = np.where(
                (folds == inner_fold) & outer_train_mask
            )[0]
            inner_train_idx = np.where(
                outer_train_mask & (folds != inner_fold)
            )[0]

            inner_model = fit_base(
                [mutation_sets[i] for i in inner_train_idx],
                labels[inner_train_idx],
            )

            for i in inner_hold_idx:
                _, b3 = score_base(mutation_sets[i], inner_model)
                inner_residual3[i] = labels[i] - b3

        if np.isnan(inner_residual3[outer_train_idx]).any():
            raise RuntimeError(
                f"Missing inner residuals for outer fold {outer_fold}."
            )

        outer_triplet_memory = build_subset_memory(
            [mutation_sets[i] for i in outer_train_idx],
            inner_residual3[outer_train_idx],
            order=3,
            shrink_lambda=2.0,
            min_support=2,
        )

        outer_base = fit_base(
            [mutation_sets[i] for i in outer_train_idx],
            labels[outer_train_idx],
        )

        for i in outer_hold_idx:
            _, b3 = score_base(mutation_sets[i], outer_base)
            triplet_delta, _, _ = adaptive_order_contribution(
                mutation_sets[i],
                outer_triplet_memory,
                3,
            )
            oof_b4[i] = b3 + triplet_delta

    if np.isnan(oof_b4).any():
        raise RuntimeError("Missing B4 OOF predictions.")

    residual4_oof = labels - oof_b4

    final_quartet_memory = build_subset_memory(
        mutation_sets,
        residual4_oof,
        order=4,
        shrink_lambda=4.0,
        min_support=2,
    )

    final_base = fit_base(mutation_sets, labels)

    return {
        "base": final_base,
        "triplet": final_triplet_memory,
        "quartet": final_quartet_memory,
        "folds": folds,
        "oof_b2": oof_b2,
        "oof_b3": oof_b3,
        "oof_b4": oof_b4,
        "residual3_oof": residual3_oof,
        "residual4_oof": residual4_oof,
    }


def score_hierarchy(mutation_set, model):
    b2, b3 = score_base(mutation_set, model["base"])

    triplet_delta, triplet_supported, triplet_confidence = (
        adaptive_order_contribution(
            mutation_set,
            model["triplet"],
            3,
        )
    )
    b4 = b3 + triplet_delta

    quartet_delta, quartet_supported, quartet_confidence = (
        adaptive_order_contribution(
            mutation_set,
            model["quartet"],
            4,
        )
    )
    b5 = b4 + quartet_delta

    return {
        "B2_ADDITIVE": float(b2),
        "B3_RAW_PAIR": float(b3),
        "B4_CROSSFIT_TRIPLET": float(b4),
        "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER": float(b5),
        "triplet_delta": float(triplet_delta),
        "quartet_delta": float(quartet_delta),
        "triplet_supported": int(triplet_supported),
        "quartet_supported": int(quartet_supported),
        "triplet_confidence": float(triplet_confidence),
        "quartet_confidence": float(quartet_confidence),
    }


def evaluate(frame, score_column, truth_column):
    truth = frame[truth_column].to_numpy(dtype=float)
    ranked = frame.sort_values(
        [score_column, "candidate_id"],
        ascending=[False, True],
    )
    top1_truth = float(ranked.iloc[0][truth_column])
    top5_truth = ranked.head(5)[truth_column].to_numpy(dtype=float)
    top50_truth = ranked.head(50)[truth_column].to_numpy(dtype=float)
    top1_cutoff = float(np.quantile(truth, 0.99))
    best = float(np.max(truth))
    worst = float(np.min(truth))

    rho = spearmanr(
        frame[score_column],
        frame[truth_column],
    ).statistic

    return {
        "spearman": float(rho),
        "top1_percentile": percentile(truth, top1_truth),
        "top5_mean_percentile": float(
            np.mean([percentile(truth, x) for x in top5_truth])
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
            [left, "candidate_id"],
            ascending=[False, True],
        ).head(k)["candidate_id"]
    )
    right_ids = set(
        frame.sort_values(
            [right, "candidate_id"],
            ascending=[False, True],
        ).head(k)["candidate_id"]
    )
    return int(len(left_ids & right_ids))


def permute_within_mutation_count(labels, mutation_sets, seed):
    labels = np.asarray(labels, dtype=float)
    counts = np.array(
        [len(mutation_set) for mutation_set in mutation_sets],
        dtype=int,
    )
    result = labels.copy()
    rng = np.random.default_rng(seed)

    for mutation_count in sorted(set(counts.tolist())):
        idx = np.where(counts == mutation_count)[0]
        result[idx] = rng.permutation(labels[idx])

    return result


def summarize_memory(memory):
    supports = np.array(
        [entry["support"] for entry in memory.values()],
        dtype=float,
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


def bootstrap_spearman_delta(
    frame,
    left_score,
    right_score,
    truth_column,
    reps=1000,
    seed=SEED,
):
    rng = np.random.default_rng(seed)
    left = frame[left_score].to_numpy(dtype=float)
    right = frame[right_score].to_numpy(dtype=float)
    truth = frame[truth_column].to_numpy(dtype=float)
    n = len(frame)

    deltas = []
    for _ in range(reps):
        idx = rng.integers(0, n, size=n)
        left_rho = spearmanr(left[idx], truth[idx]).statistic
        right_rho = spearmanr(right[idx], truth[idx]).statistic
        if np.isfinite(left_rho) and np.isfinite(right_rho):
            deltas.append(float(left_rho - right_rho))

    values = np.asarray(deltas, dtype=float)
    return {
        "reps_valid": int(len(values)),
        "mean_delta": float(np.mean(values)),
        "ci95_low": float(np.quantile(values, 0.025)),
        "ci95_high": float(np.quantile(values, 0.975)),
        "fraction_delta_gt_0": float(np.mean(values > 0)),
    }


def model_diagnostics(model):
    return {
        "main": summarize_memory(model["base"]["main"]),
        "pair": summarize_memory(model["base"]["pair"]),
        "triplet_crossfit": summarize_memory(model["triplet"]),
        "quartet_crossfit": summarize_memory(model["quartet"]),
        "oof": {
            "b3_spearman_visible": float(
                spearmanr(model["oof_b3"], model["oof_b3"] + model["residual3_oof"]).statistic
            ),
            "b4_spearman_visible": float(
                spearmanr(model["oof_b4"], model["oof_b4"] + model["residual4_oof"]).statistic
            ),
            "residual3_std": float(np.std(model["residual3_oof"])),
            "residual4_std": float(np.std(model["residual4_oof"])),
        },
        "fold_counts": {
            str(fold): int(np.sum(model["folds"] == fold))
            for fold in range(N_FOLDS)
        },
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
    data["candidate_id"] = data["mutation_set"].map(
        lambda x: ":".join(x)
    )
    data["bucket"] = data["candidate_id"].map(
        lambda x: fnv1a32(x) % 10
    )
    data["mutation_count"] = data["mutation_set"].map(len)

    visible = data[data["bucket"] <= 6].copy()
    hidden = data[data["bucket"] >= 7].copy()

    visible_sets = visible["mutation_set"].tolist()
    visible_labels = visible[truth_column].to_numpy(dtype=float)
    visible_ids = visible["candidate_id"].tolist()

    print(
        f"[V8.1] rows={len(data)} visible={len(visible)} hidden={len(hidden)}"
    )
    print("[V8.1] fitting real cross-fitted hierarchy...")

    model = fit_crossfitted_hierarchy(
        visible_sets,
        visible_labels,
        visible_ids,
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
        score_hierarchy(mutation_set, model)
        for mutation_set in eligible["mutation_set"]
    ]
    scored_frame = pd.DataFrame(scored, index=eligible.index)
    for column in scored_frame.columns:
        eligible[column] = scored_frame[column]

    print("[V8.1] fitting global permutation control...")
    global_perm_labels = np.random.default_rng(SEED).permutation(
        visible_labels
    )
    global_perm_model = fit_crossfitted_hierarchy(
        visible_sets,
        global_perm_labels,
        visible_ids,
    )
    eligible["GLOBAL_PERMUTED_B5"] = [
        score_hierarchy(mutation_set, global_perm_model)[
            "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER"
        ]
        for mutation_set in eligible["mutation_set"]
    ]

    print("[V8.1] fitting mutation-count-stratified permutation control...")
    count_perm_labels = permute_within_mutation_count(
        visible_labels,
        visible_sets,
        SEED,
    )
    count_perm_model = fit_crossfitted_hierarchy(
        visible_sets,
        count_perm_labels,
        visible_ids,
    )
    eligible["COUNT_STRATIFIED_PERMUTED_B5"] = [
        score_hierarchy(mutation_set, count_perm_model)[
            "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER"
        ]
        for mutation_set in eligible["mutation_set"]
    ]

    arms = [
        "B2_ADDITIVE",
        "B3_RAW_PAIR",
        "B4_CROSSFIT_TRIPLET",
        "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER",
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

    b3 = metrics["B3_RAW_PAIR"]
    b4 = metrics["B4_CROSSFIT_TRIPLET"]
    b5 = metrics["B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER"]

    bootstrap = bootstrap_spearman_delta(
        eligible,
        "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER",
        "B3_RAW_PAIR",
        truth_column,
        reps=1000,
    )

    result = {
        "version": "NABU_V8_1_CROSSFIT_HIGHER_ORDER_DEV",
        "status": "DEVELOPMENT_DIAGNOSTIC_NOT_FROZEN",
        "dataset": Path(csv_path).name,
        "rows": int(len(data)),
        "visible": int(len(visible)),
        "hidden": int(len(hidden)),
        "eligible_hidden": int(len(eligible)),
        "metrics": metrics,
        "deltas_vs_raw_pair": {
            "B4_minus_B3_spearman": float(
                b4["spearman"] - b3["spearman"]
            ),
            "B5_minus_B3_spearman": float(
                b5["spearman"] - b3["spearman"]
            ),
            "B5_minus_B3_top5_mean_percentile": float(
                b5["top5_mean_percentile"]
                - b3["top5_mean_percentile"]
            ),
            "B5_minus_B3_top50_enrichment": float(
                b5["top50_top1pct_enrichment"]
                - b3["top50_top1pct_enrichment"]
            ),
            "B5_minus_B3_top1_percentile": float(
                b5["top1_percentile"]
                - b3["top1_percentile"]
            ),
            "B5_minus_B3_normalized_regret": float(
                b5["normalized_top1_regret"]
                - b3["normalized_top1_regret"]
            ),
        },
        "candidate_overlap_B3_vs_B5": {
            "top1": overlap_count(
                eligible,
                "B3_RAW_PAIR",
                "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER",
                1,
            ),
            "top5": overlap_count(
                eligible,
                "B3_RAW_PAIR",
                "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER",
                5,
            ),
            "top10": overlap_count(
                eligible,
                "B3_RAW_PAIR",
                "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER",
                10,
            ),
            "top20": overlap_count(
                eligible,
                "B3_RAW_PAIR",
                "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER",
                20,
            ),
            "top50": overlap_count(
                eligible,
                "B3_RAW_PAIR",
                "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER",
                50,
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
            "mean_triplet_delta": float(
                np.mean(eligible["triplet_delta"])
            ),
            "mean_quartet_delta": float(
                np.mean(eligible["quartet_delta"])
            ),
        },
        "bootstrap_B5_minus_B3_spearman": bootstrap,
        "model_diagnostics": model_diagnostics(model),
        "by_mutation_count": by_mutation_count,
        "interpretation_boundary": (
            "PHOT development diagnostic after prior reveal. "
            "Cross-fitting tests robustness of the higher-order mechanism "
            "but does not constitute a fresh blind confirmation."
        ),
    }

    (out / "V8_1_CROSSFIT_DEV_RESULTS.json").write_text(
        json.dumps(result, indent=2),
        encoding="utf-8",
    )

    diagnostics = {
        "version": "NABU_V8_1_CROSSFIT_DIAGNOSTICS",
        "model": model_diagnostics(model),
        "higher_order_hidden_coverage": result[
            "higher_order_hidden_coverage"
        ],
        "bootstrap_B5_minus_B3_spearman": bootstrap,
    }
    (out / "V8_1_CROSSFIT_DIAGNOSTICS.json").write_text(
        json.dumps(diagnostics, indent=2),
        encoding="utf-8",
    )

    eligible.drop(columns=["mutation_set"]).to_csv(
        out / "V8_1_CROSSFIT_REVEALED_PREDICTIONS.csv",
        index=False,
    )

    top_rows = []
    for arm in arms[:4]:
        ranked = eligible.sort_values(
            [arm, "candidate_id"],
            ascending=[False, True],
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
        out / "V8_1_CROSSFIT_TOP50.csv",
        index=False,
    )

    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("csv")
    parser.add_argument(
        "--out",
        default="nabu_v8_1_crossfit_dev_results",
    )
    args = parser.parse_args()
    main(args.csv, args.out)
