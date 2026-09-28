from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from collections import Counter
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.preprocessing import PowerTransformer

HERE = Path(__file__).resolve().parent
DEV = HERE.parent / "v9_pair_transfer_dev"
sys.path.insert(0, str(DEV))

from nabu_protein.higher_order import additive_score, pair_score
from nabu_protein.v83 import NabuV83Model
from run_v9_pair_transfer_ired import derive_reference, evaluate, mutation_set

SEED = 161
EPS = 1e-12

DATASETS = {
    "TrpB": {
        "hash_type": "md5",
        "expected_hash": "a611408d2db06907a023ddc8a1d94c12",
        "expected_counts": {"fit": 8633, "validation": 2158, "test": 217507},
    },
    "IRED": {
        "hash_type": "sha256",
        "expected_hash": "aa45a2f85fb1af87b6b0e86397b3f8f292a061dc574be2536673b5bd490b0e74",
        "expected_counts": {"fit": 3746, "validation": 662, "test": 4178},
    },
}


def file_hash(path: Path, algorithm: str) -> str:
    h = hashlib.new(algorithm)
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_file(path: Path) -> str:
    return file_hash(path, "sha256")


def validation_mask(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series.dtype):
        return series.fillna(False).astype(bool)
    text = series.fillna("").astype(str).str.strip().str.lower()
    return text.isin({"true", "1", "yes", "y", "t"})


def split_frame(source: Path, dataset: str):
    cfg = DATASETS[dataset]
    observed_hash = file_hash(source, cfg["hash_type"])
    if observed_hash != cfg["expected_hash"]:
        raise RuntimeError(
            f"{dataset} hash mismatch: expected {cfg['expected_hash']}, got {observed_hash}"
        )

    frame = pd.read_csv(source, compression="gzip")
    required = {"sequence", "target", "set", "validation"}
    missing = required - set(frame.columns)
    if missing:
        raise RuntimeError(f"{dataset} missing columns: {sorted(missing)}")

    frame = frame.copy()
    frame["sequence"] = frame["sequence"].astype(str).str.strip().str.upper()
    frame["target"] = pd.to_numeric(frame["target"], errors="raise")
    if not np.isfinite(frame["target"].to_numpy(dtype=float)).all():
        raise RuntimeError(f"{dataset} has non-finite targets.")

    set_text = frame["set"].astype(str).str.strip().str.lower()
    val = validation_mask(frame["validation"])
    fit = frame[set_text.eq("train") & ~val].copy()
    validation = frame[set_text.eq("train") & val].copy()
    test = frame[set_text.eq("test")].copy()

    counts = {
        "fit": int(len(fit)),
        "validation": int(len(validation)),
        "test": int(len(test)),
    }
    if counts != cfg["expected_counts"]:
        raise RuntimeError(
            f"{dataset} split count mismatch: expected {cfg['expected_counts']}, got {counts}"
        )

    return fit, validation, test, {
        "hash_type": cfg["hash_type"],
        "source_hash": observed_hash,
        "counts": counts,
    }


def reconstruct(reference: str, mutations: tuple[str, ...]) -> str:
    chars = list(reference)
    seen = set()
    for token in mutations:
        source = token[0]
        target = token[-1]
        position = int(token[1:-1]) - 1
        if position in seen:
            raise RuntimeError(f"Duplicate position in mutation set: {mutations}")
        seen.add(position)
        if chars[position] != source:
            raise RuntimeError(
                f"Mutation source mismatch at {position + 1}: "
                f"reference={chars[position]} token={token}"
            )
        chars[position] = target
    return "".join(chars)


def unique_lookup(mutation_sets, targets, label):
    lookup = {}
    duplicates = []
    for muts, target in zip(mutation_sets, targets):
        key = tuple(muts)
        if key in lookup:
            duplicates.append(key)
        else:
            lookup[key] = float(target)
    return lookup, duplicates


def exact_components(lookup):
    if () not in lookup:
        raise RuntimeError("WT missing from fit lookup.")
    wt = float(lookup[()])
    e1 = {
        muts[0]: float(value - wt)
        for muts, value in lookup.items()
        if len(muts) == 1
    }
    e2 = {}
    for muts, value in lookup.items():
        if len(muts) != 2:
            continue
        left, right = muts
        if left in e1 and right in e1:
            e2[muts] = float(value - wt - e1[left] - e1[right])
    return wt, e1, e2


def o1_o2_predictions(mutation_sets, wt, e1, e2, shuffled_e2=None):
    n = len(mutation_sets)
    o1 = np.full(n, np.nan, dtype=float)
    o2 = np.full(n, np.nan, dtype=float)
    shuffled = np.full(n, np.nan, dtype=float)
    o1_support = np.zeros(n, dtype=bool)
    o2_support = np.zeros(n, dtype=bool)
    direct_support = np.zeros(n, dtype=bool)

    for i, muts in enumerate(mutation_sets):
        muts = tuple(muts)
        if all(m in e1 for m in muts):
            o1_support[i] = True
            o1[i] = float(wt + sum(e1[m] for m in muts))
            pairs = [tuple(pair) for pair in combinations(muts, 2)]
            if all(pair in e2 for pair in pairs):
                o2_support[i] = True
                o2[i] = float(o1[i] + sum(e2[pair] for pair in pairs))
                if shuffled_e2 is not None:
                    shuffled[i] = float(
                        o1[i] + sum(shuffled_e2[pair] for pair in pairs)
                    )

        singleton_keys = [(m,) for m in muts]
        pair_keys = [tuple(pair) for pair in combinations(muts, 2)]
        # Independent direct exact-subset support check is filled later using
        # lookup-derived keys via a caller-provided closure.
        direct_support[i] = bool(o2_support[i])

    return o1, o2, shuffled, o1_support, o2_support, direct_support


def finite_spearman(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) < 2 or np.std(a) <= EPS or np.std(b) <= EPS:
        return None
    x = float(spearmanr(a, b).statistic)
    return x if np.isfinite(x) else None


def finite_pearson(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) < 2 or np.std(a) <= EPS or np.std(b) <= EPS:
        return None
    x = float(pearsonr(a, b).statistic)
    return x if np.isfinite(x) else None


def extended_metrics(target, prediction):
    target = np.asarray(target, dtype=float)
    prediction = np.asarray(prediction, dtype=float)
    base = evaluate(target, prediction)
    err = prediction - target
    sd = float(np.std(target))
    rmse = float(np.sqrt(np.mean(err ** 2)))
    base.update(
        {
            "pearson": finite_pearson(target, prediction),
            "rmse": rmse,
            "mae": float(np.mean(np.abs(err))),
            "target_std": sd,
            "normalized_rmse": None if sd <= EPS else float(rmse / sd),
        }
    )
    return base


def masked_metrics(target, prediction, mask):
    target = np.asarray(target, dtype=float)
    prediction = np.asarray(prediction, dtype=float)
    mask = np.asarray(mask, dtype=bool)
    mask &= np.isfinite(target) & np.isfinite(prediction)
    if int(mask.sum()) < 2:
        return {"count": int(mask.sum()), "spearman": None}
    return extended_metrics(target[mask], prediction[mask])


def learned_low_order_alignment(fit_sets, fit_target, fit_ids, e1, e2):
    model = NabuV83Model().fit(
        mutation_sets=fit_sets,
        labels=fit_target,
        candidate_ids=fit_ids,
    )
    base = model.model["base"]

    shared_single = sorted(set(e1) & set(base["main"]))
    exact_single = np.asarray([e1[k] for k in shared_single], dtype=float)
    learned_single = np.asarray(
        [float(base["main"][k]["effect"]) for k in shared_single],
        dtype=float,
    )

    shared_pair = sorted(set(e2) & set(base["pair"]))
    exact_pair = np.asarray([e2[k] for k in shared_pair], dtype=float)
    learned_pair = np.asarray(
        [
            float(base["pair"][k]["mean"] * base["pair"][k]["confidence"])
            for k in shared_pair
        ],
        dtype=float,
    )

    def alignment(exact, learned):
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

    return model, {
        "single": alignment(exact_single, learned_single),
        "pair": alignment(exact_pair, learned_pair),
    }


def b2_b3_predictions(model, mutation_sets):
    base = model.model["base"]
    global_mean = float(base["global_mean"])
    b2 = np.empty(len(mutation_sets), dtype=float)
    b3 = np.empty(len(mutation_sets), dtype=float)
    main_supported = np.zeros(len(mutation_sets), dtype=bool)
    for i, muts in enumerate(mutation_sets):
        ok = all(m in base["main"] for m in muts)
        main_supported[i] = ok
        if ok:
            b2[i] = float(additive_score(muts, global_mean, base["main"]))
            b3[i] = float(pair_score(muts, b2[i], base["pair"]))
        else:
            b2[i] = global_mean
            b3[i] = global_mean
    return b2, b3, main_supported


def invariants_for_dataset(
    dataset,
    fit,
    test,
    reference,
    fit_sets,
    test_sets,
    fit_target,
    lookup,
    duplicates,
    wt,
    e1,
    e2,
):
    reconstructed = [
        reconstruct(reference, muts) == seq
        for seq, muts in zip(fit["sequence"].astype(str), fit_sets)
    ]

    fit_o1, fit_o2, _, fit_o1_support, fit_o2_support, _ = o1_o2_predictions(
        fit_sets, wt, e1, e2
    )

    order = np.asarray([len(x) for x in fit_sets], dtype=int)
    low_mask = order <= 1
    clean_double = (order == 2) & fit_o2_support

    o1_err = np.abs(fit_o1[low_mask] - fit_target[low_mask])
    o2_err = np.abs(fit_o2[clean_double] - fit_target[clean_double])

    test_o1, test_o2, _, _, test_o2_support, _ = o1_o2_predictions(
        test_sets, wt, e1, e2
    )

    # Independent support check from exact subset presence in the fit lookup.
    direct = []
    for muts in test_sets:
        singleton_ok = all((m,) in lookup for m in muts)
        pair_ok = all(tuple(pair) in lookup for pair in combinations(muts, 2))
        direct.append(singleton_ok and pair_ok)
    direct = np.asarray(direct, dtype=bool)

    synthetic = np.asarray([3.0, -1.0, 2.0, 7.0, 0.5], dtype=float)
    metric_checks = {
        "spearman_self_error": abs(float(spearmanr(synthetic, synthetic).statistic) - 1.0),
        "pearson_self_error": abs(float(pearsonr(synthetic, synthetic).statistic) - 1.0),
        "rmse_self": float(np.sqrt(np.mean((synthetic - synthetic) ** 2))),
    }

    checks = {
        "no_duplicate_fit_mutation_sets": len(duplicates) == 0,
        "reference_parses_to_empty": mutation_set(reference, reference) == (),
        "fit_sequence_roundtrip_all": bool(all(reconstructed)),
        "o1_fit_reconstruction_max_abs_le_1e12": bool(
            len(o1_err) > 0 and float(np.max(o1_err)) <= 1e-12
        ),
        "o2_clean_double_reconstruction_max_abs_le_1e12": bool(
            len(o2_err) > 0 and float(np.max(o2_err)) <= 1e-12
        ),
        "o2_support_bookkeeping_exact": bool(np.array_equal(test_o2_support, direct)),
        "spearman_self_within_1e12": metric_checks["spearman_self_error"] <= 1e-12,
        "pearson_self_within_1e12": metric_checks["pearson_self_error"] <= 1e-12,
        "rmse_self_within_1e12": metric_checks["rmse_self"] <= 1e-12,
    }

    return {
        "dataset": dataset,
        "checks": checks,
        "all_required_pass": bool(all(checks.values())),
        "fit_duplicate_mutation_set_count": int(len(duplicates)),
        "fit_roundtrip_fail_count": int(len(reconstructed) - sum(reconstructed)),
        "o1_fit_reconstruction_count": int(np.sum(low_mask)),
        "o1_fit_reconstruction_max_abs_error": (
            None if len(o1_err) == 0 else float(np.max(o1_err))
        ),
        "o2_clean_double_reconstruction_count": int(np.sum(clean_double)),
        "o2_clean_double_reconstruction_max_abs_error": (
            None if len(o2_err) == 0 else float(np.max(o2_err))
        ),
        "test_o2_support_count": int(np.sum(test_o2_support)),
        "direct_support_count": int(np.sum(direct)),
        "support_disagreement_count": int(np.sum(test_o2_support != direct)),
        "metric_self_tests": metric_checks,
    }


class IdentityTransform:
    name = "identity"

    def fit(self, y):
        return self

    def transform(self, y):
        return np.asarray(y, dtype=float)


class YeoJohnsonTransform:
    name = "yeo_johnson"

    def __init__(self):
        self.pt = PowerTransformer(method="yeo-johnson", standardize=False)

    def fit(self, y):
        self.pt.fit(np.asarray(y, dtype=float).reshape(-1, 1))
        return self

    def transform(self, y):
        return self.pt.transform(
            np.asarray(y, dtype=float).reshape(-1, 1)
        ).reshape(-1)


class RobustAsinhTransform:
    name = "robust_asinh"

    def fit(self, y):
        y = np.asarray(y, dtype=float)
        self.center = float(np.median(y))
        mad = float(np.median(np.abs(y - self.center)))
        self.scale_source = "mad"
        if mad <= EPS:
            mad = float(np.std(y))
            self.scale_source = "std_fallback"
        if mad <= EPS:
            mad = 1.0
            self.scale_source = "unit_fallback"
        self.scale = mad
        return self

    def transform(self, y):
        y = np.asarray(y, dtype=float)
        return np.arcsinh((y - self.center) / self.scale)


def transformed_o2_diagnostics(
    fit_sets,
    fit_target,
    test_sets,
    test_target,
    test_orders,
):
    output = {}
    for transform in (IdentityTransform(), YeoJohnsonTransform(), RobustAsinhTransform()):
        transform.fit(fit_target)
        fit_y = transform.transform(fit_target)
        test_y = transform.transform(test_target)
        lookup, duplicates = unique_lookup(fit_sets, fit_y, transform.name)
        if duplicates:
            raise RuntimeError(
                f"Duplicate transformed fit mutation sets for {transform.name}"
            )
        wt, e1, e2 = exact_components(lookup)

        keys = sorted(e2)
        values = np.asarray([e2[k] for k in keys], dtype=float)
        shuffled_values = np.random.default_rng(SEED).permutation(values)
        shuffled_e2 = {k: float(v) for k, v in zip(keys, shuffled_values)}

        o1, o2, shuffled, _, o2_support, _ = o1_o2_predictions(
            test_sets, wt, e1, e2, shuffled_e2
        )
        per_order = {}
        for order in sorted(set(test_orders.tolist())):
            mask = (test_orders == order) & o2_support
            if int(np.sum(mask)) >= 2:
                per_order[str(order)] = {
                    "support_count": int(np.sum(mask)),
                    "O1": masked_metrics(test_y, o1, mask),
                    "O2": masked_metrics(test_y, o2, mask),
                    "SHUFFLED_O2": masked_metrics(test_y, shuffled, mask),
                }
        meta = {}
        if isinstance(transform, YeoJohnsonTransform):
            meta["lambda"] = float(transform.pt.lambdas_[0])
        if isinstance(transform, RobustAsinhTransform):
            meta.update(
                {
                    "center": float(transform.center),
                    "scale": float(transform.scale),
                    "scale_source": transform.scale_source,
                }
            )
        output[transform.name] = {"fit_only_parameters": meta, "by_order": per_order}
    return output


def pair_context_diagnostics(test_sets, test_target, test_orders, wt, e1, e2):
    triple_e3 = {}
    base_pairs = []
    conditional_pairs = []
    supported_triples = 0

    for muts, target, order in zip(test_sets, test_target, test_orders):
        if order != 3:
            continue
        if not all(m in e1 for m in muts):
            continue
        pairs = [tuple(pair) for pair in combinations(muts, 2)]
        if not all(pair in e2 for pair in pairs):
            continue
        o2 = float(
            wt
            + sum(e1[m] for m in muts)
            + sum(e2[pair] for pair in pairs)
        )
        e3 = float(target - o2)
        triple_e3[tuple(muts)] = e3
        supported_triples += 1
        for pair in pairs:
            base_pairs.append(float(e2[pair]))
            conditional_pairs.append(float(e2[pair] + e3))

    base_pairs = np.asarray(base_pairs, dtype=float)
    conditional_pairs = np.asarray(conditional_pairs, dtype=float)
    eligible = (np.abs(base_pairs) > EPS) & (np.abs(conditional_pairs) > EPS)
    sign_flip = (
        None
        if int(np.sum(eligible)) == 0
        else float(
            np.mean(
                np.sign(base_pairs[eligible])
                != np.sign(conditional_pairs[eligible])
            )
        )
    )

    e3_values = np.asarray(list(triple_e3.values()), dtype=float)
    triple_targets = np.asarray(
        [
            target
            for muts, target, order in zip(test_sets, test_target, test_orders)
            if order == 3 and tuple(muts) in triple_e3
        ],
        dtype=float,
    )
    target_sd = float(np.std(triple_targets)) if len(triple_targets) else 0.0
    e3_rms = (
        None
        if len(e3_values) == 0
        else float(np.sqrt(np.mean(e3_values ** 2)))
    )

    return {
        "supported_triple_count": int(supported_triples),
        "pair_context_instance_count": int(len(base_pairs)),
        "sign_flip_rate": sign_flip,
        "spearman_base_vs_conditional": finite_spearman(
            base_pairs, conditional_pairs
        ),
        "pearson_base_vs_conditional": finite_pearson(
            base_pairs, conditional_pairs
        ),
        "e3_rms": e3_rms,
        "triple_target_std": target_sd,
        "e3_rms_over_target_std": (
            None
            if e3_rms is None or target_sd <= EPS
            else float(e3_rms / target_sd)
        ),
    }


def dataset_diagnostics(dataset, source):
    fit, validation, test, source_info = split_frame(source, dataset)

    fit_sequences = fit["sequence"].astype(str).tolist()
    test_sequences = test["sequence"].astype(str).tolist()
    fit_target = fit["target"].to_numpy(dtype=float)
    test_target = test["target"].to_numpy(dtype=float)

    reference = derive_reference(fit_sequences)
    fit_sets = [mutation_set(seq, reference) for seq in fit_sequences]
    test_sets = [mutation_set(seq, reference) for seq in test_sequences]
    fit_orders = np.asarray([len(x) for x in fit_sets], dtype=int)
    test_orders = np.asarray([len(x) for x in test_sets], dtype=int)

    lookup, duplicates = unique_lookup(fit_sets, fit_target, f"{dataset}_fit")
    wt, e1, e2 = exact_components(lookup)

    invariants = invariants_for_dataset(
        dataset,
        fit,
        test,
        reference,
        fit_sets,
        test_sets,
        fit_target,
        lookup,
        duplicates,
        wt,
        e1,
        e2,
    )

    model, representation = learned_low_order_alignment(
        fit_sets,
        fit_target,
        fit["sequence"].astype(str).tolist(),
        e1,
        e2,
    )
    b2, b3, b2_support = b2_b3_predictions(model, test_sets)

    keys = sorted(e2)
    pair_values = np.asarray([e2[k] for k in keys], dtype=float)
    shuffled_values = np.random.default_rng(SEED).permutation(pair_values)
    shuffled_e2 = {k: float(v) for k, v in zip(keys, shuffled_values)}

    o1, o2, shuffled_o2, _, o2_support, _ = o1_o2_predictions(
        test_sets, wt, e1, e2, shuffled_e2
    )

    by_order = {}
    coverage_bias = {}
    for order in sorted(set(test_orders.tolist())):
        order_mask = test_orders == order
        supported = order_mask & o2_support
        unsupported = order_mask & ~o2_support

        block = {
            "count": int(np.sum(order_mask)),
            "o2_supported_count": int(np.sum(supported)),
            "o2_unsupported_count": int(np.sum(unsupported)),
        }
        if int(np.sum(supported)) >= 2:
            block["O1_supported"] = masked_metrics(test_target, o1, supported)
            block["O2"] = masked_metrics(test_target, o2, supported)
            block["SHUFFLED_O2"] = masked_metrics(
                test_target, shuffled_o2, supported
            )
            block["B2_supported"] = masked_metrics(
                test_target, b2, supported
            )
            block["B3_supported"] = masked_metrics(
                test_target, b3, supported
            )
        if int(np.sum(unsupported)) >= 2:
            block["B2_unsupported"] = masked_metrics(
                test_target, b2, unsupported
            )

        if int(np.sum(supported)) >= 100 and int(np.sum(unsupported)) >= 100:
            s = block["B2_supported"]["spearman"]
            u = block["B2_unsupported"]["spearman"]
            coverage_bias[str(order)] = {
                "evaluable": True,
                "b2_supported_spearman": s,
                "b2_unsupported_spearman": u,
                "absolute_difference": (
                    None if s is None or u is None else float(abs(s - u))
                ),
                "target_mean_supported": float(np.mean(test_target[supported])),
                "target_mean_unsupported": float(np.mean(test_target[unsupported])),
                "target_std_supported": float(np.std(test_target[supported])),
                "target_std_unsupported": float(np.std(test_target[unsupported])),
            }
        else:
            coverage_bias[str(order)] = {
                "evaluable": False,
                "supported_count": int(np.sum(supported)),
                "unsupported_count": int(np.sum(unsupported)),
            }

        by_order[str(order)] = block

    transforms = transformed_o2_diagnostics(
        fit_sets, fit_target, test_sets, test_target, test_orders
    )
    context = pair_context_diagnostics(
        test_sets, test_target, test_orders, wt, e1, e2
    )

    return {
        "dataset": dataset,
        "source": source_info,
        "reference_length": int(len(reference)),
        "reference_sha256": hashlib.sha256(reference.encode("utf-8")).hexdigest(),
        "fit_order_distribution": {
            str(k): int(v) for k, v in sorted(Counter(fit_orders).items())
        },
        "test_order_distribution": {
            str(k): int(v) for k, v in sorted(Counter(test_orders).items())
        },
        "invariants": invariants,
        "representation": representation,
        "coverage_bias": coverage_bias,
        "context_dependence": context,
        "by_order": by_order,
        "target_scale_transforms": transforms,
        "b2_main_supported_count": int(np.sum(b2_support)),
        "exact_component_counts": {
            "single": int(len(e1)),
            "clean_pair": int(len(e2)),
        },
    }


def classify_matrix(diag):
    trpb = diag["TrpB"]
    ired = diag["IRED"]

    # H0 implementation/data
    both_invariants = (
        trpb["invariants"]["all_required_pass"]
        and ired["invariants"]["all_required_pass"]
    )
    h0 = {
        "status": "REJECTED" if both_invariants else "SUPPORTED",
        "reason": (
            "All required decomposition/data/metric invariants passed on both datasets."
            if both_invariants
            else "At least one required invariant failed."
        ),
    }

    # H1 representation
    def rep_pass(d):
        s = d["representation"]["single"]
        p = d["representation"]["pair"]
        return (
            s["count"] >= 30
            and p["count"] >= 30
            and s["spearman"] is not None
            and p["spearman"] is not None
            and s["spearman"] >= 0.50
            and p["spearman"] >= 0.50
        )

    tr_rep = rep_pass(trpb)
    ir_rep = rep_pass(ired)
    if tr_rep and ir_rep:
        h1_status = "REJECTED"
    elif (not tr_rep) and (not ir_rep):
        h1_status = "SUPPORTED"
    else:
        h1_status = "UNRESOLVED"
    h1 = {
        "status": h1_status,
        "trpb_pass": tr_rep,
        "ired_pass": ir_rep,
        "details": {
            "TrpB": trpb["representation"],
            "IRED": ired["representation"],
        },
    }

    # H2 coverage bias
    common_orders = sorted(
        set(trpb["coverage_bias"]) & set(ired["coverage_bias"]),
        key=int,
    )
    comparable = []
    for order in common_orders:
        a = trpb["coverage_bias"][order]
        b = ired["coverage_bias"][order]
        if a.get("evaluable") and b.get("evaluable"):
            comparable.append(
                (
                    order,
                    a["absolute_difference"],
                    b["absolute_difference"],
                    a["b2_supported_spearman"] - a["b2_unsupported_spearman"],
                    b["b2_supported_spearman"] - b["b2_unsupported_spearman"],
                )
            )
    if comparable:
        same_dir_material = any(
            aa is not None
            and bb is not None
            and aa >= 0.10
            and bb >= 0.10
            and np.sign(da) == np.sign(db)
            for _, aa, bb, da, db in comparable
        )
        both_small = all(
            aa is not None and bb is not None and aa < 0.10 and bb < 0.10
            for _, aa, bb, _, _ in comparable
        )
        h2_status = (
            "SUPPORTED"
            if same_dir_material
            else ("REJECTED" if both_small else "UNRESOLVED")
        )
    else:
        h2_status = "UNRESOLVED"
    h2 = {
        "status": h2_status,
        "matched_evaluable_orders": [
            {
                "order": order,
                "trpb_abs_difference": aa,
                "ired_abs_difference": bb,
                "trpb_signed_difference": da,
                "ired_signed_difference": db,
            }
            for order, aa, bb, da, db in comparable
        ],
    }

    # H3 context-free pair portability
    def portability_evidence(d, order):
        block = d["by_order"].get(order, {})
        if block.get("o2_supported_count", 0) < 100:
            return None
        for key in ("O1_supported", "O2", "SHUFFLED_O2"):
            if key not in block or block[key].get("spearman") is None:
                return None
        o1s = block["O1_supported"]["spearman"]
        o2s = block["O2"]["spearman"]
        shs = block["SHUFFLED_O2"]["spearman"]
        return {
            "o1": o1s,
            "o2": o2s,
            "shuffled": shs,
            "o2_minus_o1": float(o2s - o1s),
            "o2_minus_shuffled": float(o2s - shs),
            "non_portability": bool(
                (o2s - o1s < 0.05) or (o2s <= shs + 0.02)
            ),
            "portable": bool(
                (o2s - o1s >= 0.05) and (o2s > shs + 0.02)
            ),
        }

    common_test_orders = sorted(
        set(trpb["by_order"]) & set(ired["by_order"]), key=int
    )
    pair_common = []
    for order in common_test_orders:
        a = portability_evidence(trpb, order)
        b = portability_evidence(ired, order)
        if a is not None and b is not None:
            pair_common.append((order, a, b))
    if any(a["non_portability"] and b["non_portability"] for _, a, b in pair_common):
        h3_status = "SUPPORTED"
    elif pair_common and all(a["portable"] and b["portable"] for _, a, b in pair_common):
        h3_status = "REJECTED"
    else:
        h3_status = "UNRESOLVED"
    h3 = {
        "status": h3_status,
        "matched_evaluable_orders": [
            {"order": order, "TrpB": a, "IRED": b}
            for order, a, b in pair_common
        ],
    }

    # H4 context dependence
    def context_material(d):
        c = d["context_dependence"]
        if c["supported_triple_count"] < 100:
            return None
        if c["sign_flip_rate"] is None or c["e3_rms_over_target_std"] is None:
            return None
        return bool(
            c["sign_flip_rate"] >= 0.20
            and c["e3_rms_over_target_std"] >= 0.50
        )

    tc = context_material(trpb)
    ic = context_material(ired)
    if tc is True and ic is True:
        h4_status = "SUPPORTED"
    elif tc is False and ic is False:
        h4_status = "REJECTED"
    else:
        h4_status = "UNRESOLVED"
    h4 = {
        "status": h4_status,
        "TrpB": trpb["context_dependence"],
        "IRED": ired["context_dependence"],
    }

    # H5 target-scale nonlinearity
    def transform_improvements(d, order):
        identity = (
            d["target_scale_transforms"]["identity"]["by_order"].get(order)
        )
        if identity is None or identity["O2"]["spearman"] is None:
            return {}
        base_s = identity["O2"]["spearman"]
        base_nrmse = identity["O2"]["normalized_rmse"]
        out = {}
        for name in ("yeo_johnson", "robust_asinh"):
            block = d["target_scale_transforms"][name]["by_order"].get(order)
            if block is None or block["O2"]["spearman"] is None:
                continue
            s_gain = float(block["O2"]["spearman"] - base_s)
            if base_nrmse is None or base_nrmse <= EPS:
                rmse_reduction = None
            else:
                rmse_reduction = float(
                    (base_nrmse - block["O2"]["normalized_rmse"]) / base_nrmse
                )
            out[name] = {
                "spearman_gain": s_gain,
                "normalized_rmse_reduction_fraction": rmse_reduction,
                "material": bool(
                    s_gain >= 0.10
                    and rmse_reduction is not None
                    and rmse_reduction >= 0.20
                ),
                "reaches_either": bool(
                    s_gain >= 0.10
                    or (
                        rmse_reduction is not None
                        and rmse_reduction >= 0.20
                    )
                ),
            }
        return out

    transform_common = []
    for order in common_test_orders:
        a = transform_improvements(trpb, order)
        b = transform_improvements(ired, order)
        if a and b:
            transform_common.append((order, a, b))

    supported_transform = False
    for order, a, b in transform_common:
        for name in ("yeo_johnson", "robust_asinh"):
            if (
                name in a
                and name in b
                and a[name]["material"]
                and b[name]["material"]
            ):
                supported_transform = True

    any_either = any(
        item["reaches_either"]
        for _, a, b in transform_common
        for item in list(a.values()) + list(b.values())
    )
    if supported_transform:
        h5_status = "SUPPORTED"
    elif transform_common and not any_either:
        h5_status = "REJECTED"
    else:
        h5_status = "UNRESOLVED"
    h5 = {
        "status": h5_status,
        "matched_evaluable_orders": [
            {"order": order, "TrpB": a, "IRED": b}
            for order, a, b in transform_common
        ],
    }

    h6 = {
        "status": "UNRESOLVED",
        "reason": (
            "Stage-1 does not claim fundamental non-identifiability. "
            "A separate strict-CV higher-order residual structure test is required."
        ),
    }

    hypotheses = {
        "H0_IMPLEMENTATION_OR_DATA_SEMANTICS": h0,
        "H1_REPRESENTATION_OR_ESTIMATION_FAILURE": h1,
        "H2_COVERAGE_OR_SUPPORT_BIAS": h2,
        "H3_CONTEXT_FREE_PAIR_NON_PORTABILITY": h3,
        "H4_CONTEXT_DEPENDENT_HIGHER_ORDER_INTERACTION": h4,
        "H5_OBSERVATION_SCALE_OR_COMPOSITION_NONLINEARITY": h5,
        "H6_LOW_ORDER_INFORMATION_INSUFFICIENCY": h6,
    }

    supported = [
        key for key, value in hypotheses.items() if value["status"] == "SUPPORTED"
    ]
    rejected = [
        key for key, value in hypotheses.items() if value["status"] == "REJECTED"
    ]

    # Conservative global rule: Stage-1 cannot declare root cause while H6 is
    # intentionally unresolved, and any cross-dataset evidence gap blocks it.
    root_status = "ROOT_CAUSE_NOT_YET_IDENTIFIED"

    return {
        "version": "NABU_ROOT_CAUSE_FALSIFICATION_STAGE1_V1",
        "root_cause_status": root_status,
        "hypotheses": hypotheses,
        "supported_hypotheses": supported,
        "rejected_hypotheses": rejected,
        "unresolved_hypotheses": [
            key
            for key, value in hypotheses.items()
            if value["status"] == "UNRESOLVED"
        ],
        "nucb_consumed": False,
        "phase3_opened": False,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--trpb", required=True)
    parser.add_argument("--ired", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    trpb_path = Path(args.trpb)
    ired_path = Path(args.ired)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    diagnostics = {
        "TrpB": dataset_diagnostics("TrpB", trpb_path),
        "IRED": dataset_diagnostics("IRED", ired_path),
    }

    matrix = classify_matrix(diagnostics)

    (out / "TRPB_INVARIANTS.json").write_text(
        json.dumps(diagnostics["TrpB"]["invariants"], indent=2),
        encoding="utf-8",
    )
    (out / "IRED_INVARIANTS.json").write_text(
        json.dumps(diagnostics["IRED"]["invariants"], indent=2),
        encoding="utf-8",
    )
    (out / "DATASET_DIAGNOSTICS.json").write_text(
        json.dumps(diagnostics, indent=2),
        encoding="utf-8",
    )
    (out / "FALSIFICATION_MATRIX.json").write_text(
        json.dumps(matrix, indent=2),
        encoding="utf-8",
    )

    manifest = {
        "version": "NABU_ROOT_CAUSE_FALSIFICATION_RUN_MANIFEST_V1",
        "seed": SEED,
        "inputs": {
            "TrpB": {
                "path_name": trpb_path.name,
                "sha256": sha256_file(trpb_path),
            },
            "IRED": {
                "path_name": ired_path.name,
                "sha256": sha256_file(ired_path),
            },
        },
        "nucb_consumed": False,
        "phase3_opened": False,
    }
    manifest_path = out / "RUN_MANIFEST.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    hashes = {}
    for name in (
        "TRPB_INVARIANTS.json",
        "IRED_INVARIANTS.json",
        "DATASET_DIAGNOSTICS.json",
        "FALSIFICATION_MATRIX.json",
        "RUN_MANIFEST.json",
    ):
        hashes[name] = sha256_file(out / name)
    (out / "OUTPUT_HASHES.json").write_text(
        json.dumps(hashes, indent=2),
        encoding="utf-8",
    )

    print(json.dumps(matrix, indent=2))


if __name__ == "__main__":
    main()
