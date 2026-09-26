from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

HERE = Path(__file__).resolve().parent
DEV = HERE.parent / "v9_pair_transfer_dev"
sys.path.insert(0, str(DEV))

from nabu_protein.higher_order import additive_score, pair_score
from nabu_protein.v83 import NabuV83Model
from run_v9_pair_transfer_ired import derive_reference, evaluate, mutation_set

SEED = 161
MATERIAL_RANK_GAIN = 0.10
NEGLIGIBLE_PAIR_IDENTITY_GAIN = 0.02
MATERIAL_RESIDUAL = 0.50
MATERIAL_SIGN_INSTABILITY = 0.20
REPRESENTATION_GAP = 0.10
MATERIAL_ERROR_ALIGNMENT = 0.50
STRONG_ORACLE_SPEARMAN = 0.50
EPS = 1e-12

PRIMARY_PRIORITY = [
    "MISSING_THIRD_ORDER_INFORMATION_IS_CAUSAL_FOR_QUADS",
    "FOURTH_ORDER_INTERACTION_IS_MATERIAL",
    "CONTEXT_DEPENDENT_HIGHER_ORDER_INTERACTION",
    "PAIR_REPRESENTATION_OR_ESTIMATION_FAILURE",
    "CONTEXT_FREE_PAIR_IDENTITY_NOT_SUFFICIENT",
    "B2_FAILURE_TRACKS_MISSING_HIGHER_ORDER_TERMS",
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def finite_spearman(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) < 2 or np.std(a) <= EPS or np.std(b) <= EPS:
        return None
    value = float(spearmanr(a, b).statistic)
    return value if np.isfinite(value) else None


def finite_pearson(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) < 2 or np.std(a) <= EPS or np.std(b) <= EPS:
        return None
    value = float(pearsonr(a, b).statistic)
    return value if np.isfinite(value) else None


def extended_metrics(target, prediction):
    target = np.asarray(target, dtype=float)
    prediction = np.asarray(prediction, dtype=float)
    base = evaluate(target, prediction)
    error = prediction - target
    base.update(
        {
            "pearson": finite_pearson(target, prediction),
            "rmse": float(np.sqrt(np.mean(error ** 2))),
            "mae": float(np.mean(np.abs(error))),
            "target_mean": float(np.mean(target)),
            "target_std": float(np.std(target)),
            "prediction_mean": float(np.mean(prediction)),
            "prediction_std": float(np.std(prediction)),
        }
    )
    return base


def masked_metrics(target, prediction, mask):
    target = np.asarray(target, dtype=float)
    prediction = np.asarray(prediction, dtype=float)
    mask = np.asarray(mask, dtype=bool)
    mask &= np.isfinite(target) & np.isfinite(prediction)
    if int(mask.sum()) < 2:
        return {
            "count": int(mask.sum()),
            "spearman": None,
            "ndcg": None,
            "top1_percent_hits": None,
            "rmse": None,
            "mae": None,
            "pearson": None,
        }
    return extended_metrics(target[mask], prediction[mask])


def distribution_stats(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return {"count": 0}
    return {
        "count": int(len(values)),
        "mean": float(np.mean(values)),
        "std": float(np.std(values)),
        "rms": float(np.sqrt(np.mean(values ** 2))),
        "mean_abs": float(np.mean(np.abs(values))),
        "q05": float(np.quantile(values, 0.05)),
        "q25": float(np.quantile(values, 0.25)),
        "median": float(np.quantile(values, 0.50)),
        "q75": float(np.quantile(values, 0.75)),
        "q95": float(np.quantile(values, 0.95)),
    }


def normalized_residual_stats(residual, target):
    stats = distribution_stats(residual)
    target = np.asarray(target, dtype=float)
    target = target[np.isfinite(target)]
    target_std = float(np.std(target)) if len(target) else 0.0
    stats["target_std"] = target_std
    stats["rms_over_target_std"] = (
        None
        if target_std <= EPS or stats.get("count", 0) == 0
        else float(stats["rms"] / target_std)
    )
    return stats


def unique_target_lookup(mutation_sets, targets, name):
    lookup = {}
    for muts, target in zip(mutation_sets, targets):
        key = tuple(muts)
        if key in lookup:
            raise RuntimeError(f"Duplicate mutation set in {name}: {key}")
        lookup[key] = float(target)
    return lookup


def pair_context_stability(triple_sets, triple_e3, e2):
    base = []
    conditional = []
    for muts in triple_sets:
        e3 = triple_e3.get(tuple(muts))
        if e3 is None:
            continue
        for pair in combinations(muts, 2):
            pair = tuple(pair)
            if pair not in e2:
                raise RuntimeError("Internal support inconsistency in pair context audit.")
            base.append(float(e2[pair]))
            conditional.append(float(e2[pair] + e3))

    base = np.asarray(base, dtype=float)
    conditional = np.asarray(conditional, dtype=float)
    nonzero = (np.abs(base) > EPS) & (np.abs(conditional) > EPS)
    if int(nonzero.sum()) == 0:
        sign_flip = None
    else:
        sign_flip = float(
            np.mean(np.sign(base[nonzero]) != np.sign(conditional[nonzero]))
        )

    return {
        "pair_context_instances": int(len(base)),
        "sign_flip_eligible_count": int(nonzero.sum()),
        "sign_flip_rate": sign_flip,
        "spearman_base_vs_conditional": finite_spearman(base, conditional),
        "pearson_base_vs_conditional": finite_pearson(base, conditional),
        "base_stats": distribution_stats(base),
        "conditional_stats": distribution_stats(conditional),
    }


def alignment_report(exact, learned):
    exact = np.asarray(exact, dtype=float)
    learned = np.asarray(learned, dtype=float)
    return {
        "count": int(len(exact)),
        "spearman": finite_spearman(exact, learned),
        "pearson": finite_pearson(exact, learned),
        "exact_std": float(np.std(exact)) if len(exact) else None,
        "learned_std": float(np.std(learned)) if len(learned) else None,
        "learned_to_exact_std_ratio": (
            None
            if len(exact) == 0 or np.std(exact) <= EPS
            else float(np.std(learned) / np.std(exact))
        ),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fit_csv_gz")
    parser.add_argument("test_masked_csv_gz")
    parser.add_argument("test_reveal_csv_gz")
    parser.add_argument("--source-manifest", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    fit_path = Path(args.fit_csv_gz)
    masked_path = Path(args.test_masked_csv_gz)
    reveal_path = Path(args.test_reveal_csv_gz)
    manifest_path = Path(args.source_manifest)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    fit = pd.read_csv(fit_path, compression="gzip")
    masked = pd.read_csv(masked_path, compression="gzip")
    reveal = pd.read_csv(reveal_path, compression="gzip")
    source_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    test = masked.merge(reveal, on="row_id", how="inner", validate="one_to_one")
    if len(test) != len(masked) or len(test) != len(reveal):
        raise RuntimeError("Masked/reveal merge mismatch.")

    fit_sequences = fit["sequence"].astype(str).str.strip().str.upper().tolist()
    test_sequences = test["sequence"].astype(str).str.strip().str.upper().tolist()
    fit_target = fit["target"].to_numpy(dtype=float)
    test_target = test["target"].to_numpy(dtype=float)
    fit_ids = fit["row_id"].astype(str).tolist()

    reference = derive_reference(fit_sequences)
    fit_sets = [mutation_set(sequence, reference) for sequence in fit_sequences]
    test_sets = [mutation_set(sequence, reference) for sequence in test_sequences]

    fit_orders = np.asarray([len(x) for x in fit_sets], dtype=int)
    test_orders = np.asarray([len(x) for x in test_sets], dtype=int)
    if set(fit_orders.tolist()) - {0, 1, 2}:
        raise RuntimeError("Fit contains mutation order above 2.")
    if set(test_orders.tolist()) - {3, 4}:
        raise RuntimeError("Test contains mutation order outside 3/4.")

    fit_lookup = unique_target_lookup(fit_sets, fit_target, "fit")
    _ = unique_target_lookup(test_sets, test_target, "test")

    if () not in fit_lookup:
        raise RuntimeError("No exact WT row in fit.")
    wt = float(fit_lookup[()])

    singleton_y = {
        muts[0]: value
        for muts, value in fit_lookup.items()
        if len(muts) == 1
    }
    e1 = {mutation: float(value - wt) for mutation, value in singleton_y.items()}

    e2 = {}
    raw_double_count = 0
    excluded_double_missing_single = 0
    for muts, value in fit_lookup.items():
        if len(muts) != 2:
            continue
        raw_double_count += 1
        left, right = muts
        if left not in e1 or right not in e1:
            excluded_double_missing_single += 1
            continue
        e2[muts] = float(value - wt - e1[left] - e1[right])

    pair_keys = sorted(e2)
    pair_values = np.asarray([e2[key] for key in pair_keys], dtype=float)
    shuffled_values = np.random.default_rng(SEED).permutation(pair_values)
    e2_shuffled = {
        key: float(value)
        for key, value in zip(pair_keys, shuffled_values)
    }

    model = NabuV83Model().fit(
        mutation_sets=fit_sets,
        labels=fit_target,
        candidate_ids=fit_ids,
    )
    base = model.model["base"]
    global_mean = float(base["global_mean"])

    n = len(test)
    b2 = np.empty(n, dtype=float)
    b3 = np.empty(n, dtype=float)
    b2_supported = np.zeros(n, dtype=bool)
    b3_exact_pair_count = np.zeros(n, dtype=int)

    o1 = np.full(n, np.nan, dtype=float)
    o2 = np.full(n, np.nan, dtype=float)
    o2_shuffled = np.full(n, np.nan, dtype=float)
    o1_supported = np.zeros(n, dtype=bool)
    o2_supported = np.zeros(n, dtype=bool)

    for i, muts in enumerate(test_sets):
        main_ok = all(mutation in base["main"] for mutation in muts)
        b2_supported[i] = main_ok
        if main_ok:
            b2[i] = float(additive_score(muts, global_mean, base["main"]))
            b3[i] = float(pair_score(muts, b2[i], base["pair"]))
            b3_exact_pair_count[i] = int(
                sum(tuple(pair) in base["pair"] for pair in combinations(muts, 2))
            )
        else:
            b2[i] = global_mean
            b3[i] = global_mean
            b3_exact_pair_count[i] = 0

        if all(mutation in e1 for mutation in muts):
            o1_supported[i] = True
            o1[i] = float(wt + sum(e1[mutation] for mutation in muts))

            pairs = [tuple(pair) for pair in combinations(muts, 2)]
            if all(pair in e2 for pair in pairs):
                o2_supported[i] = True
                o2[i] = float(o1[i] + sum(e2[pair] for pair in pairs))
                o2_shuffled[i] = float(
                    o1[i] + sum(e2_shuffled[pair] for pair in pairs)
                )

    e3 = np.full(n, np.nan, dtype=float)
    triple_e3 = {}
    for i, muts in enumerate(test_sets):
        if test_orders[i] == 3 and o2_supported[i]:
            value = float(test_target[i] - o2[i])
            e3[i] = value
            triple_e3[tuple(muts)] = value

    o3 = np.full(n, np.nan, dtype=float)
    o3_supported = np.zeros(n, dtype=bool)
    e4 = np.full(n, np.nan, dtype=float)
    for i, muts in enumerate(test_sets):
        if test_orders[i] != 4 or not o2_supported[i]:
            continue
        triples = [tuple(x) for x in combinations(muts, 3)]
        if all(triple in triple_e3 for triple in triples):
            o3_supported[i] = True
            o3[i] = float(o2[i] + sum(triple_e3[triple] for triple in triples))
            e4[i] = float(test_target[i] - o3[i])

    rows = pd.DataFrame(
        {
            "row_id": test["row_id"].astype(str),
            "mutation_count": test_orders,
            "target": test_target,
            "b2_supported": b2_supported,
            "b3_exact_pair_count": b3_exact_pair_count,
            "B2": b2,
            "B3_RAW_PAIR": b3,
            "o1_supported": o1_supported,
            "EXACT_O1": o1,
            "o2_supported": o2_supported,
            "EXACT_O2": o2,
            "SHUFFLED_O2": o2_shuffled,
            "O2_RESIDUAL": test_target - o2,
            "E3": e3,
            "o3_supported": o3_supported,
            "EXACT_O3": o3,
            "E4": e4,
        }
    )

    metrics = {}
    pair_identity_gain = {}
    representation_gap = {}

    for order in (3, 4):
        order_mask = test_orders == order
        metrics[str(order)] = {
            "all_rows": {
                "B2": masked_metrics(test_target, b2, order_mask),
                "B3_RAW_PAIR": masked_metrics(test_target, b3, order_mask),
            },
            "o1_supported": {
                "support_count": int(np.sum(order_mask & o1_supported)),
                "support_fraction": float(
                    np.mean(o1_supported[order_mask])
                ),
                "EXACT_O1": masked_metrics(
                    test_target, o1, order_mask & o1_supported
                ),
            },
            "o2_supported": {
                "support_count": int(np.sum(order_mask & o2_supported)),
                "support_fraction": float(
                    np.mean(o2_supported[order_mask])
                ),
                "EXACT_O2": masked_metrics(
                    test_target, o2, order_mask & o2_supported
                ),
                "SHUFFLED_O2": masked_metrics(
                    test_target, o2_shuffled, order_mask & o2_supported
                ),
                "B2_same_rows": masked_metrics(
                    test_target, b2, order_mask & o2_supported
                ),
                "B3_same_rows": masked_metrics(
                    test_target, b3, order_mask & o2_supported
                ),
            },
        }

        exact_s = metrics[str(order)]["o2_supported"]["EXACT_O2"]["spearman"]
        shuffled_s = metrics[str(order)]["o2_supported"]["SHUFFLED_O2"]["spearman"]
        b3_s = metrics[str(order)]["o2_supported"]["B3_same_rows"]["spearman"]
        pair_identity_gain[str(order)] = (
            None
            if exact_s is None or shuffled_s is None
            else float(exact_s - shuffled_s)
        )
        representation_gap[str(order)] = (
            None
            if exact_s is None or b3_s is None
            else float(exact_s - b3_s)
        )

    quad_o3_mask = (test_orders == 4) & o3_supported
    metrics["4"]["o3_supported"] = {
        "support_count": int(np.sum(quad_o3_mask)),
        "support_fraction_of_quads": float(np.mean(o3_supported[test_orders == 4])),
        "EXACT_O2_same_rows": masked_metrics(test_target, o2, quad_o3_mask),
        "EXACT_O3": masked_metrics(test_target, o3, quad_o3_mask),
        "B2_same_rows": masked_metrics(test_target, b2, quad_o3_mask),
        "B3_same_rows": masked_metrics(test_target, b3, quad_o3_mask),
    }

    triple_o2_mask = (test_orders == 3) & o2_supported
    quad_o2_mask = (test_orders == 4) & o2_supported

    interaction_magnitude = {
        "e1": distribution_stats(list(e1.values())),
        "e2": distribution_stats(list(e2.values())),
        "e3_triples": normalized_residual_stats(
            e3[triple_o2_mask], test_target[triple_o2_mask]
        ),
        "o2_residual_quads": normalized_residual_stats(
            (test_target - o2)[quad_o2_mask],
            test_target[quad_o2_mask],
        ),
        "e4_quads": normalized_residual_stats(
            e4[quad_o3_mask], test_target[quad_o3_mask]
        ),
    }

    stability = pair_context_stability(
        [test_sets[i] for i in np.where(triple_o2_mask)[0]],
        triple_e3,
        e2,
    )

    shared_mutations = sorted(set(e1) & set(base["main"]))
    exact_main = [e1[key] for key in shared_mutations]
    learned_main = [float(base["main"][key]["effect"]) for key in shared_mutations]

    shared_pairs = sorted(set(e2) & set(base["pair"]))
    exact_pair = [e2[key] for key in shared_pairs]
    learned_pair = [
        float(base["pair"][key]["mean"] * base["pair"][key]["confidence"])
        for key in shared_pairs
    ]

    representation = {
        "b2_main_vs_exact_single": alignment_report(exact_main, learned_main),
        "b3_pair_contribution_vs_exact_e2": alignment_report(
            exact_pair, learned_pair
        ),
        "b2_main_memory_count": int(len(base["main"])),
        "b3_pair_memory_count": int(len(base["pair"])),
        "exact_single_count": int(len(e1)),
        "exact_e2_count": int(len(e2)),
        "raw_double_count": int(raw_double_count),
        "excluded_double_missing_single": int(excluded_double_missing_single),
    }

    error_attribution = {}
    alignment_values = []
    for order, mask in [(3, triple_o2_mask), (4, quad_o2_mask)]:
        residual = (test_target - o2)[mask]
        b2_error = (test_target - b2)[mask]
        b3_error = (test_target - b3)[mask]
        block = {
            "count": int(np.sum(mask)),
            "o2_residual_vs_b2_error_spearman": finite_spearman(
                residual, b2_error
            ),
            "o2_residual_vs_b3_error_spearman": finite_spearman(
                residual, b3_error
            ),
            "o2_residual_vs_b2_error_pearson": finite_pearson(
                residual, b2_error
            ),
            "o2_residual_vs_b3_error_pearson": finite_pearson(
                residual, b3_error
            ),
        }
        error_attribution[str(order)] = block
        for key in (
            "o2_residual_vs_b2_error_spearman",
            "o2_residual_vs_b3_error_spearman",
        ):
            if block[key] is not None:
                alignment_values.append(abs(float(block[key])))

    if int(np.sum(quad_o3_mask)) >= 2:
        e4_values = e4[quad_o3_mask]
        b2_error = (test_target - b2)[quad_o3_mask]
        b3_error = (test_target - b3)[quad_o3_mask]
        quad_e4_block = {
            "count": int(np.sum(quad_o3_mask)),
            "e4_vs_b2_error_spearman": finite_spearman(e4_values, b2_error),
            "e4_vs_b3_error_spearman": finite_spearman(e4_values, b3_error),
            "e4_vs_b2_error_pearson": finite_pearson(e4_values, b2_error),
            "e4_vs_b3_error_pearson": finite_pearson(e4_values, b3_error),
        }
    else:
        quad_e4_block = {"count": int(np.sum(quad_o3_mask))}
    error_attribution["quad_e4"] = quad_e4_block
    for key in ("e4_vs_b2_error_spearman", "e4_vs_b3_error_spearman"):
        value = quad_e4_block.get(key)
        if value is not None:
            alignment_values.append(abs(float(value)))

    triggered = []

    representation_failure = False
    for order in ("3", "4"):
        exact_s = metrics[order]["o2_supported"]["EXACT_O2"]["spearman"]
        gap = representation_gap[order]
        if (
            exact_s is not None
            and gap is not None
            and exact_s >= STRONG_ORACLE_SPEARMAN
            and gap >= REPRESENTATION_GAP
        ):
            representation_failure = True
    if representation_failure:
        triggered.append("PAIR_REPRESENTATION_OR_ESTIMATION_FAILURE")

    e3_ratio = interaction_magnitude["e3_triples"].get("rms_over_target_std")
    sign_flip = stability.get("sign_flip_rate")
    if (
        e3_ratio is not None
        and e3_ratio >= MATERIAL_RESIDUAL
        and sign_flip is not None
        and sign_flip >= MATERIAL_SIGN_INSTABILITY
    ):
        triggered.append("CONTEXT_DEPENDENT_HIGHER_ORDER_INTERACTION")

    o3_s = metrics["4"]["o3_supported"]["EXACT_O3"]["spearman"]
    o2_same_s = metrics["4"]["o3_supported"]["EXACT_O2_same_rows"]["spearman"]
    o3_gain = (
        None
        if o3_s is None or o2_same_s is None
        else float(o3_s - o2_same_s)
    )
    if o3_gain is not None and o3_gain >= MATERIAL_RANK_GAIN:
        triggered.append("MISSING_THIRD_ORDER_INFORMATION_IS_CAUSAL_FOR_QUADS")

    e4_ratio = interaction_magnitude["e4_quads"].get("rms_over_target_std")
    if e4_ratio is not None and e4_ratio >= MATERIAL_RESIDUAL:
        triggered.append("FOURTH_ORDER_INTERACTION_IS_MATERIAL")

    valid_pair_gains = [
        value for value in pair_identity_gain.values() if value is not None
    ]
    if valid_pair_gains and all(
        value <= NEGLIGIBLE_PAIR_IDENTITY_GAIN for value in valid_pair_gains
    ):
        triggered.append("CONTEXT_FREE_PAIR_IDENTITY_NOT_SUFFICIENT")

    max_error_alignment = max(alignment_values) if alignment_values else None
    if (
        max_error_alignment is not None
        and max_error_alignment >= MATERIAL_ERROR_ALIGNMENT
    ):
        triggered.append("B2_FAILURE_TRACKS_MISSING_HIGHER_ORDER_TERMS")

    triggered = list(dict.fromkeys(triggered))
    primary = next(
        (cause for cause in PRIMARY_PRIORITY if cause in triggered),
        "ROOT_CAUSE_NOT_RESOLVED_BY_THIS_AUDIT",
    )

    decision = {
        "primary_cause": primary,
        "triggered_causes": triggered,
        "o3_minus_o2_spearman_same_supported_quads": o3_gain,
        "pair_identity_spearman_gain_over_shuffle": pair_identity_gain,
        "oracle_minus_b3_spearman_gap": representation_gap,
        "max_abs_error_alignment_spearman": max_error_alignment,
        "thresholds": {
            "material_rank_gain": MATERIAL_RANK_GAIN,
            "negligible_pair_identity_gain": NEGLIGIBLE_PAIR_IDENTITY_GAIN,
            "material_residual_rms_over_target_std": MATERIAL_RESIDUAL,
            "material_sign_instability": MATERIAL_SIGN_INSTABILITY,
            "representation_gap": REPRESENTATION_GAP,
            "material_error_alignment": MATERIAL_ERROR_ALIGNMENT,
            "strong_oracle_spearman": STRONG_ORACLE_SPEARMAN,
        },
    }

    report = {
        "version": "NABU_PHASE2D_ROOT_CAUSE_AUDIT_V1",
        "status": "REVEALED_TRPB_DEVELOPMENT_DIAGNOSTIC",
        "scientific_claim": False,
        "phase3_opened": False,
        "nucb_consumed": False,
        "dataset": "FLIP2_TrpB_two_to_many",
        "counts": {
            "fit": int(len(fit)),
            "test": int(len(test)),
            "fit_order_distribution": {
                str(order): int(np.sum(fit_orders == order))
                for order in sorted(set(fit_orders.tolist()))
            },
            "test_order_distribution": {
                str(order): int(np.sum(test_orders == order))
                for order in sorted(set(test_orders.tolist()))
            },
        },
        "coverage": {
            "o1_supported_total": int(np.sum(o1_supported)),
            "o2_supported_total": int(np.sum(o2_supported)),
            "o3_supported_quads": int(np.sum(o3_supported)),
            "b2_main_supported_total": int(np.sum(b2_supported)),
        },
        "metrics": metrics,
        "interaction_magnitude": interaction_magnitude,
        "pair_context_stability": stability,
        "representation": representation,
        "error_attribution": error_attribution,
        "decision": decision,
        "interpretation_boundary": {
            "context_free_o2_failure_means": (
                "WT+single+pair context-free composition is insufficient on the tested region."
            ),
            "does_not_prove": (
                "It does not prove that every possible model using low-order data must fail."
            ),
            "identifiability_note": (
                "If measured e3 materially recovers quad ranking, those exact third-order "
                "coefficients are not directly observed in WT/single/double labels and require "
                "an additional prior, representation, or additional measurements."
            ),
        },
    }

    components = {
        "version": "NABU_INTERACTION_COMPONENT_SUMMARY_V1",
        "wt": wt,
        "counts": {
            "e1": int(len(e1)),
            "e2": int(len(e2)),
            "e3_supported_triples": int(len(triple_e3)),
            "e4_supported_quads": int(np.sum(np.isfinite(e4))),
        },
        "magnitude": interaction_magnitude,
        "pair_context_stability": stability,
        "representation": representation,
    }

    run_manifest = {
        "version": "NABU_ROOT_CAUSE_RUN_MANIFEST_V1",
        "seed": SEED,
        "source_manifest": source_manifest,
        "input_sha256": {
            "fit": sha256_file(fit_path),
            "test_masked": sha256_file(masked_path),
            "test_reveal": sha256_file(reveal_path),
            "source_manifest": sha256_file(manifest_path),
        },
        "reference_sha256": hashlib.sha256(reference.encode("utf-8")).hexdigest(),
        "reference_length": int(len(reference)),
        "decision_thresholds": decision["thresholds"],
        "nucb_consumed": False,
        "phase3_opened": False,
    }

    report_path = out / "ROOT_CAUSE_REPORT.json"
    rows_path = out / "ROW_LEVEL_AUDIT.csv.gz"
    components_path = out / "INTERACTION_COMPONENTS.json"
    manifest_out_path = out / "RUN_MANIFEST.json"

    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    rows.to_csv(\n        rows_path,\n        index=False,\n        compression={"method": "gzip", "mtime": 0},\n    )
    components_path.write_text(json.dumps(components, indent=2), encoding="utf-8")
    manifest_out_path.write_text(
        json.dumps(run_manifest, indent=2), encoding="utf-8"
    )

    output_hashes = {
        "ROOT_CAUSE_REPORT.json": sha256_file(report_path),
        "ROW_LEVEL_AUDIT.csv.gz": sha256_file(rows_path),
        "INTERACTION_COMPONENTS.json": sha256_file(components_path),
        "RUN_MANIFEST.json": sha256_file(manifest_out_path),
    }
    (out / "OUTPUT_HASHES.json").write_text(
        json.dumps(output_hashes, indent=2), encoding="utf-8"
    )

    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
