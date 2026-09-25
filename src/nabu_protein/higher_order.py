"""
Canonical NABU V8.3 higher-order hierarchy.

This module packages the frozen Phase-1 B2/B3/B4/B5 implementation without
changing its mathematics. Source of truth:
freeze/nabu-v8-3-dual-objective-router-validated
commit c8afdcd7698231d95a39ef13a3fe22b6e3f507c5
"""

from __future__ import annotations

import math
import re
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

SEED = 161
N_FOLDS = 5
FROZEN_CORE_COMMIT = "c8afdcd7698231d95a39ef13a3fe22b6e3f507c5"


def fnv1a32(text: str) -> int:
    h = 2166136261
    for byte in text.encode("utf-8"):
        h = ((h ^ byte) * 16777619) & 0xFFFFFFFF
    return h


def parse_mutations(value):
    tokens = re.findall(r"[A-Z][0-9]+[A-Z*]", str(value))
    return tuple(sorted(tokens, key=lambda token: (int(token[1:-1]), token)))


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


def model_diagnostics(model):
    return {
        "main": summarize_memory(model["base"]["main"]),
        "pair": summarize_memory(model["base"]["pair"]),
        "triplet_crossfit": summarize_memory(model["triplet"]),
        "quartet_crossfit": summarize_memory(model["quartet"]),
        "oof": {
            "b3_spearman_visible": float(
                spearmanr(
                    model["oof_b3"],
                    model["oof_b3"] + model["residual3_oof"],
                ).statistic
            ),
            "b4_spearman_visible": float(
                spearmanr(
                    model["oof_b4"],
                    model["oof_b4"] + model["residual4_oof"],
                ).statistic
            ),
            "residual3_std": float(np.std(model["residual3_oof"])),
            "residual4_std": float(np.std(model["residual4_oof"])),
        },
        "fold_counts": {
            str(fold): int(np.sum(model["folds"] == fold))
            for fold in range(N_FOLDS)
        },
    }
