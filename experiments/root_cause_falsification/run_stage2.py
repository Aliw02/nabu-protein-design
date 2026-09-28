from __future__ import annotations

import argparse
import hashlib
import json
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from run_falsification import (
    EPS,
    SEED,
    dataset_diagnostics,
    exact_components,
    finite_spearman,
    learned_low_order_alignment,
    masked_metrics,
    o1_o2_predictions,
    pair_context_diagnostics,
    transformed_o2_diagnostics,
    unique_lookup,
)
from nabu_protein.higher_order import additive_score, fnv1a32, pair_score
from nabu_protein.v83 import NabuV83Model

GB1_WT = "VDGV"
GB1_GIT_BLOB_SHA = "a1e9f5146a9441210eaf88accfa660f02d80cbb7"
GB1_EXPECTED_SIZE = 2993024
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def git_blob_sha(path: Path) -> str:
    data = path.read_bytes()
    payload = f"blob {len(data)}\0".encode("utf-8") + data
    return hashlib.sha1(payload).hexdigest()


def gb1_mutation_set(genotype: str) -> tuple[str, ...]:
    genotype = str(genotype).strip().upper()
    if len(genotype) != 4:
        raise ValueError(f"GB1 genotype length != 4: {genotype!r}")
    if not set(genotype).issubset(STANDARD_AA):
        raise ValueError(f"GB1 non-standard amino acid: {genotype!r}")
    tokens = []
    for i, (source, target) in enumerate(zip(GB1_WT, genotype), start=1):
        if source != target:
            tokens.append(f"{source}{i}{target}")
    return tuple(tokens)


def reconstruct_gb1(muts: tuple[str, ...]) -> str:
    chars = list(GB1_WT)
    for token in muts:
        pos = int(token[1:-1]) - 1
        source = token[0]
        target = token[-1]
        if chars[pos] != source:
            raise RuntimeError(f"GB1 source mismatch for {token}")
        chars[pos] = target
    return "".join(chars)


def b2_b3_predictions(model, mutation_sets):
    base = model.model["base"]
    mean = float(base["global_mean"])
    b2 = np.empty(len(mutation_sets), dtype=float)
    b3 = np.empty(len(mutation_sets), dtype=float)
    support = np.zeros(len(mutation_sets), dtype=bool)
    for i, muts in enumerate(mutation_sets):
        ok = all(m in base["main"] for m in muts)
        support[i] = ok
        if ok:
            b2[i] = float(additive_score(muts, mean, base["main"]))
            b3[i] = float(pair_score(muts, b2[i], base["pair"]))
        else:
            b2[i] = mean
            b3[i] = mean
    return b2, b3, support


def gb1_diagnostics(path: Path):
    if path.stat().st_size != GB1_EXPECTED_SIZE:
        raise RuntimeError(
            f"GB1 size drift: expected {GB1_EXPECTED_SIZE}, got {path.stat().st_size}"
        )
    actual_blob = git_blob_sha(path)
    if actual_blob != GB1_GIT_BLOB_SHA:
        raise RuntimeError(
            f"GB1 Git blob SHA mismatch: expected {GB1_GIT_BLOB_SHA}, got {actual_blob}"
        )

    frame = pd.read_csv(path)
    if set(["variant", "fitness"]) - set(frame.columns):
        raise RuntimeError("GB1 source columns missing.")

    frame = frame[["variant", "fitness"]].dropna().copy()
    frame["variant"] = frame["variant"].astype(str).str.strip().str.upper()
    frame["fitness"] = pd.to_numeric(frame["fitness"], errors="raise")
    frame = frame[
        frame["variant"].str.len().eq(4)
        & frame["variant"].map(lambda x: set(x).issubset(STANDARD_AA))
    ].copy()

    if frame["variant"].duplicated().any():
        raise RuntimeError("GB1 duplicate variants.")
    if GB1_WT not in set(frame["variant"]):
        raise RuntimeError("GB1 WT VDGV missing.")

    wt_fitness = float(frame.loc[frame["variant"].eq(GB1_WT), "fitness"].iloc[0])
    if abs(wt_fitness - 1.0) > 1e-12:
        raise RuntimeError(f"GB1 WT fitness drift: {wt_fitness}")

    frame["mutation_set"] = frame["variant"].map(gb1_mutation_set)
    frame["mutation_order"] = frame["mutation_set"].map(len)

    fit = frame[frame["mutation_order"].le(2)].copy()
    test = frame[frame["mutation_order"].ge(3)].copy()
    if len(fit) == 0 or len(test) == 0:
        raise RuntimeError("GB1 diagnostic split empty.")

    fit_sets = fit["mutation_set"].tolist()
    test_sets = test["mutation_set"].tolist()
    fit_target = fit["fitness"].to_numpy(dtype=float)
    test_target = test["fitness"].to_numpy(dtype=float)
    test_orders = test["mutation_order"].to_numpy(dtype=int)

    lookup, duplicates = unique_lookup(fit_sets, fit_target, "GB1_fit")
    if duplicates:
        raise RuntimeError("GB1 low-order mutation-set duplicates.")
    wt, e1, e2 = exact_components(lookup)

    # Invariants
    roundtrip_fail = int(
        sum(
            reconstruct_gb1(muts) != genotype
            for genotype, muts in zip(frame["variant"], frame["mutation_set"])
        )
    )
    fit_o1, fit_o2, _, fit_o1_support, fit_o2_support, _ = o1_o2_predictions(
        fit_sets, wt, e1, e2
    )
    fit_orders = fit["mutation_order"].to_numpy(dtype=int)
    low_mask = fit_orders <= 1
    clean_double = (fit_orders == 2) & fit_o2_support
    o1_err = np.abs(fit_o1[low_mask] - fit_target[low_mask])
    o2_err = np.abs(fit_o2[clean_double] - fit_target[clean_double])

    invariants = {
        "git_blob_sha_matches": actual_blob == GB1_GIT_BLOB_SHA,
        "file_size_matches": path.stat().st_size == GB1_EXPECTED_SIZE,
        "wt_present": True,
        "wt_fitness_equals_1": abs(wt_fitness - 1.0) <= 1e-12,
        "all_variant_roundtrips": roundtrip_fail == 0,
        "o1_max_abs_error_le_1e12": (
            len(o1_err) > 0 and float(np.max(o1_err)) <= 1e-12
        ),
        "o2_max_abs_error_le_1e12": (
            len(o2_err) > 0 and float(np.max(o2_err)) <= 1e-12
        ),
    }
    invariants["all_required_pass"] = bool(all(invariants.values()))

    model, representation = learned_low_order_alignment(
        fit_sets,
        fit_target,
        fit["variant"].tolist(),
        e1,
        e2,
    )
    b2, b3, b2_support = b2_b3_predictions(model, test_sets)

    keys = sorted(e2)
    values = np.asarray([e2[k] for k in keys], dtype=float)
    shuffled_values = np.random.default_rng(SEED).permutation(values)
    shuffled_e2 = {k: float(v) for k, v in zip(keys, shuffled_values)}
    o1, o2, shuffled_o2, _, o2_support, _ = o1_o2_predictions(
        test_sets, wt, e1, e2, shuffled_e2
    )

    by_order = {}
    coverage_bias = {}
    for order in sorted(set(test_orders.tolist())):
        om = test_orders == order
        sup = om & o2_support
        unsup = om & ~o2_support
        block = {
            "count": int(np.sum(om)),
            "o2_supported_count": int(np.sum(sup)),
            "o2_unsupported_count": int(np.sum(unsup)),
        }
        if int(np.sum(sup)) >= 2:
            block["O1_supported"] = masked_metrics(test_target, o1, sup)
            block["O2"] = masked_metrics(test_target, o2, sup)
            block["SHUFFLED_O2"] = masked_metrics(test_target, shuffled_o2, sup)
            block["B2_supported"] = masked_metrics(test_target, b2, sup)
            block["B3_supported"] = masked_metrics(test_target, b3, sup)
        if int(np.sum(unsup)) >= 2:
            block["B2_unsupported"] = masked_metrics(test_target, b2, unsup)

        if int(np.sum(sup)) >= 100 and int(np.sum(unsup)) >= 100:
            s = block["B2_supported"]["spearman"]
            u = block["B2_unsupported"]["spearman"]
            coverage_bias[str(order)] = {
                "evaluable": True,
                "b2_supported_spearman": s,
                "b2_unsupported_spearman": u,
                "absolute_difference": float(abs(s - u)),
                "signed_difference": float(s - u),
            }
        else:
            coverage_bias[str(order)] = {
                "evaluable": False,
                "supported_count": int(np.sum(sup)),
                "unsupported_count": int(np.sum(unsup)),
            }
        by_order[str(order)] = block

    transforms = transformed_o2_diagnostics(
        fit_sets, fit_target, test_sets, test_target, test_orders
    )
    context = pair_context_diagnostics(
        test_sets, test_target, test_orders, wt, e1, e2
    )

    return {
        "dataset": "GB1_Wu2016",
        "source": {
            "git_blob_sha": actual_blob,
            "sha256": sha256_file(path),
            "size": int(path.stat().st_size),
            "wt": GB1_WT,
        },
        "counts": {
            "all": int(len(frame)),
            "fit_low_order": int(len(fit)),
            "test_higher_order": int(len(test)),
        },
        "fit_order_distribution": {
            str(k): int(np.sum(fit_orders == k))
            for k in sorted(set(fit_orders.tolist()))
        },
        "test_order_distribution": {
            str(k): int(np.sum(test_orders == k))
            for k in sorted(set(test_orders.tolist()))
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
        "_h6": {
            "candidate_ids": test["variant"].tolist(),
            "mutation_sets": test_sets,
            "target": test_target.tolist(),
            "orders": test_orders.tolist(),
            "wt": wt,
            "e1": e1,
            "e2": {"|".join(k): v for k, v in e2.items()},
        },
    }


def build_h6_rows(candidate_ids, mutation_sets, target, orders, wt, e1, e2):
    features = []
    residual = []
    ids = []
    kept_orders = []

    for cid, muts, y, order in zip(candidate_ids, mutation_sets, target, orders):
        if not all(m in e1 for m in muts):
            continue
        pairs = [tuple(pair) for pair in combinations(muts, 2)]
        if not all(pair in e2 for pair in pairs):
            continue

        e1v = np.asarray([e1[m] for m in muts], dtype=float)
        e2v = np.asarray([e2[p] for p in pairs], dtype=float)
        o1 = float(wt + np.sum(e1v))
        o2 = float(o1 + np.sum(e2v))
        features.append(
            [
                float(order),
                o1,
                o2,
                float(np.sum(e1v)),
                float(np.mean(np.abs(e1v))),
                float(np.max(np.abs(e1v))),
                float(np.std(e1v)),
                float(np.sum(e2v)),
                float(np.mean(np.abs(e2v))),
                float(np.max(np.abs(e2v))),
                float(np.std(e2v)),
                float(np.min(e2v)),
                float(np.max(e2v)),
                float(np.mean(e2v > 0.0)),
            ]
        )
        residual.append(float(y - o2))
        ids.append(str(cid))
        kept_orders.append(int(order))

    return (
        np.asarray(features, dtype=float),
        np.asarray(residual, dtype=float),
        np.asarray(ids, dtype=str),
        np.asarray(kept_orders, dtype=int),
    )


def h6_cv(dataset_name, candidate_ids, mutation_sets, target, orders, wt, e1, e2):
    X, y, ids, kept_orders = build_h6_rows(
        candidate_ids, mutation_sets, target, orders, wt, e1, e2
    )
    if len(y) < 500:
        return {
            "dataset": dataset_name,
            "evaluable": False,
            "supported_rows": int(len(y)),
        }

    folds = np.asarray(
        [
            fnv1a32(cid + "|NABU_ROOT_CAUSE_H6|161") % 5
            for cid in ids
        ],
        dtype=int,
    )
    if set(folds.tolist()) != set(range(5)):
        raise RuntimeError(f"{dataset_name} H6 did not populate all folds.")

    rng = np.random.default_rng(SEED)
    shuffled_y = rng.permutation(y)
    pred = np.full(len(y), np.nan, dtype=float)
    shuffled_pred = np.full(len(y), np.nan, dtype=float)

    for fold in range(5):
        train = folds != fold
        hold = folds == fold

        real_model = RandomForestRegressor(
            n_estimators=200,
            min_samples_leaf=20,
            max_features="sqrt",
            random_state=SEED,
            n_jobs=1,
        )
        real_model.fit(X[train], y[train])
        pred[hold] = real_model.predict(X[hold])

        shuffled_model = RandomForestRegressor(
            n_estimators=200,
            min_samples_leaf=20,
            max_features="sqrt",
            random_state=SEED,
            n_jobs=1,
        )
        shuffled_model.fit(X[train], shuffled_y[train])
        shuffled_pred[hold] = shuffled_model.predict(X[hold])

    pooled_real = finite_spearman(y, pred)
    pooled_shuffle = finite_spearman(y, shuffled_pred)
    by_order = {}
    for order in sorted(set(kept_orders.tolist())):
        mask = kept_orders == order
        if int(np.sum(mask)) >= 100:
            by_order[str(order)] = {
                "count": int(np.sum(mask)),
                "real_oof_spearman": finite_spearman(y[mask], pred[mask]),
                "shuffled_control_spearman": finite_spearman(
                    y[mask], shuffled_pred[mask]
                ),
            }

    return {
        "dataset": dataset_name,
        "evaluable": True,
        "supported_rows": int(len(y)),
        "fold_counts": {
            str(f): int(np.sum(folds == f)) for f in range(5)
        },
        "pooled_real_oof_spearman": pooled_real,
        "pooled_shuffled_control_spearman": pooled_shuffle,
        "real_minus_shuffled": (
            None
            if pooled_real is None or pooled_shuffle is None
            else float(pooled_real - pooled_shuffle)
        ),
        "by_order": by_order,
    }


def unpack_gb1_h6(gb1):
    h = gb1["_h6"]
    e2 = {}
    for key, value in h["e2"].items():
        e2[tuple(key.split("|"))] = float(value)
    return (
        h["candidate_ids"],
        [tuple(x) for x in h["mutation_sets"]],
        np.asarray(h["target"], dtype=float),
        np.asarray(h["orders"], dtype=int),
        float(h["wt"]),
        {k: float(v) for k, v in h["e1"].items()},
        e2,
    )


def trpb_h6_inputs(path: Path):
    fit = pd.read_csv(path, compression="gzip")
    required = {"sequence", "target", "set", "validation"}
    if required - set(fit.columns):
        raise RuntimeError("TrpB source columns missing.")
    fit["sequence"] = fit["sequence"].astype(str).str.strip().str.upper()
    fit["target"] = pd.to_numeric(fit["target"], errors="raise")
    set_text = fit["set"].astype(str).str.strip().str.lower()
    val_text = fit["validation"].fillna("").astype(str).str.lower()
    val = val_text.isin({"true", "1", "yes", "y", "t"})
    low = fit[set_text.eq("train") & ~val].copy()
    high = fit[set_text.eq("test")].copy()

    from run_falsification import split_frame
    low2, _, high2, _ = split_frame(path, "TrpB")
    low = low2
    high = high2

    from run_v9_pair_transfer_ired import derive_reference, mutation_set
    reference = derive_reference(low["sequence"].astype(str).tolist())
    low_sets = [
        mutation_set(s, reference) for s in low["sequence"].astype(str)
    ]
    high_sets = [
        mutation_set(s, reference) for s in high["sequence"].astype(str)
    ]
    low_y = low["target"].to_numpy(dtype=float)
    high_y = high["target"].to_numpy(dtype=float)
    lookup, dup = unique_lookup(low_sets, low_y, "TrpB_H6")
    if dup:
        raise RuntimeError("TrpB H6 duplicate low-order mutation sets.")
    wt, e1, e2 = exact_components(lookup)
    orders = np.asarray([len(x) for x in high_sets], dtype=int)
    return (
        high["sequence"].astype(str).tolist(),
        high_sets,
        high_y,
        orders,
        wt,
        e1,
        e2,
    )


def portability_item(d, order):
    block = d["by_order"].get(order, {})
    if block.get("o2_supported_count", 0) < 100:
        return None
    needed = ("O1_supported", "O2", "SHUFFLED_O2")
    if any(key not in block for key in needed):
        return None
    vals = [block[key].get("spearman") for key in needed]
    if any(v is None for v in vals):
        return None
    o1s, o2s, shs = vals
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


def transform_improvements(d, order):
    identity = d["target_scale_transforms"]["identity"]["by_order"].get(order)
    if identity is None or identity["O2"]["spearman"] is None:
        return {}
    base_s = identity["O2"]["spearman"]
    base_rmse = identity["O2"]["normalized_rmse"]
    result = {}
    for name in ("yeo_johnson", "robust_asinh"):
        block = d["target_scale_transforms"][name]["by_order"].get(order)
        if block is None or block["O2"]["spearman"] is None:
            continue
        gain = float(block["O2"]["spearman"] - base_s)
        reduction = (
            None
            if base_rmse is None or base_rmse <= EPS
            else float(
                (base_rmse - block["O2"]["normalized_rmse"]) / base_rmse
            )
        )
        result[name] = {
            "spearman_gain": gain,
            "normalized_rmse_reduction_fraction": reduction,
            "material": bool(
                gain >= 0.10
                and reduction is not None
                and reduction >= 0.20
            ),
        }
    return result


def stage2_matrix(trpb, gb1, h6):
    # H2
    h2_matches = []
    for order in sorted(set(trpb["coverage_bias"]) & set(gb1["coverage_bias"]), key=int):
        a = trpb["coverage_bias"][order]
        b = gb1["coverage_bias"][order]
        if a.get("evaluable") and b.get("evaluable"):
            da = float(
                a["b2_supported_spearman"] - a["b2_unsupported_spearman"]
            )
            db = float(
                b["b2_supported_spearman"] - b["b2_unsupported_spearman"]
            )
            h2_matches.append(
                {
                    "order": order,
                    "TrpB_abs": a["absolute_difference"],
                    "GB1_abs": b["absolute_difference"],
                    "TrpB_signed": da,
                    "GB1_signed": db,
                }
            )
    if any(
        x["TrpB_abs"] >= 0.10
        and x["GB1_abs"] >= 0.10
        and np.sign(x["TrpB_signed"]) == np.sign(x["GB1_signed"])
        for x in h2_matches
    ):
        h2_status = "SUPPORTED"
    elif h2_matches and all(
        x["TrpB_abs"] < 0.10 and x["GB1_abs"] < 0.10 for x in h2_matches
    ):
        h2_status = "REJECTED"
    else:
        h2_status = "UNRESOLVED"

    # H3
    h3_matches = []
    for order in sorted(set(trpb["by_order"]) & set(gb1["by_order"]), key=int):
        a = portability_item(trpb, order)
        b = portability_item(gb1, order)
        if a is not None and b is not None:
            h3_matches.append({"order": order, "TrpB": a, "GB1": b})
    if any(x["TrpB"]["non_portability"] and x["GB1"]["non_portability"] for x in h3_matches):
        h3_status = "SUPPORTED"
    elif h3_matches and all(x["TrpB"]["portable"] and x["GB1"]["portable"] for x in h3_matches):
        h3_status = "REJECTED"
    else:
        h3_status = "UNRESOLVED"

    # H4
    def material_context(d):
        c = d["context_dependence"]
        if c["supported_triple_count"] < 100:
            return None
        if c["sign_flip_rate"] is None or c["e3_rms_over_target_std"] is None:
            return None
        return bool(
            c["sign_flip_rate"] >= 0.20
            and c["e3_rms_over_target_std"] >= 0.50
        )
    tc = material_context(trpb)
    gc = material_context(gb1)
    if tc is True and gc is True:
        h4_status = "SUPPORTED"
    elif tc is False and gc is False:
        h4_status = "REJECTED"
    else:
        h4_status = "UNRESOLVED"

    # H5
    h5_matches = []
    supported = False
    for order in sorted(set(trpb["by_order"]) & set(gb1["by_order"]), key=int):
        a = transform_improvements(trpb, order)
        b = transform_improvements(gb1, order)
        if a and b:
            h5_matches.append({"order": order, "TrpB": a, "GB1": b})
            for name in ("yeo_johnson", "robust_asinh"):
                if (
                    name in a
                    and name in b
                    and a[name]["material"]
                    and b[name]["material"]
                ):
                    supported = True
    if supported:
        h5_status = "SUPPORTED"
    else:
        any_material = any(
            item["material"]
            for x in h5_matches
            for side in ("TrpB", "GB1")
            for item in x[side].values()
        )
        h5_status = "UNRESOLVED" if any_material else "REJECTED"

    # H6
    t = h6["TrpB"]
    g = h6["GB1"]
    h6_reject_strong_insuff = bool(
        t.get("evaluable")
        and g.get("evaluable")
        and t.get("pooled_real_oof_spearman") is not None
        and g.get("pooled_real_oof_spearman") is not None
        and t["pooled_real_oof_spearman"] >= 0.30
        and g["pooled_real_oof_spearman"] >= 0.30
        and t["real_minus_shuffled"] >= 0.20
        and g["real_minus_shuffled"] >= 0.20
    )
    h6_status = "REJECTED" if h6_reject_strong_insuff else "UNRESOLVED"

    hypotheses = {
        "H2_COVERAGE_OR_SUPPORT_BIAS": {
            "status": h2_status,
            "matched_orders": h2_matches,
        },
        "H3_CONTEXT_FREE_PAIR_NON_PORTABILITY": {
            "status": h3_status,
            "matched_orders": h3_matches,
        },
        "H4_CONTEXT_DEPENDENT_HIGHER_ORDER_INTERACTION": {
            "status": h4_status,
            "TrpB": trpb["context_dependence"],
            "GB1": gb1["context_dependence"],
        },
        "H5_OBSERVATION_SCALE_OR_COMPOSITION_NONLINEARITY": {
            "status": h5_status,
            "matched_orders": h5_matches,
        },
        "H6_LOW_ORDER_INFORMATION_INSUFFICIENCY": {
            "status": h6_status,
            "interpretation": (
                "REJECTED means only that a strong no-usable-structure claim is rejected; "
                "it does not prove full identifiability from low-order data."
            ),
            "TrpB": t,
            "GB1": g,
        },
    }

    return {
        "version": "NABU_ROOT_CAUSE_FALSIFICATION_STAGE2_V1",
        "root_cause_status": "ROOT_CAUSE_NOT_YET_IDENTIFIED",
        "hypotheses": hypotheses,
        "supported_hypotheses": [
            k for k, v in hypotheses.items() if v["status"] == "SUPPORTED"
        ],
        "rejected_hypotheses": [
            k for k, v in hypotheses.items() if v["status"] == "REJECTED"
        ],
        "unresolved_hypotheses": [
            k for k, v in hypotheses.items() if v["status"] == "UNRESOLVED"
        ],
        "nucb_consumed": False,
        "phase3_opened": False,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--trpb", required=True)
    parser.add_argument("--gb1", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    trpb = dataset_diagnostics("TrpB", Path(args.trpb))
    gb1 = gb1_diagnostics(Path(args.gb1))

    trpb_h6 = h6_cv("TrpB", *trpb_h6_inputs(Path(args.trpb)))
    gb1_inputs = unpack_gb1_h6(gb1)
    gb1_h6 = h6_cv("GB1", *gb1_inputs)
    h6 = {"TrpB": trpb_h6, "GB1": gb1_h6}

    matrix = stage2_matrix(trpb, gb1, h6)

    # Remove internal payload before serialization.
    gb1_public = dict(gb1)
    gb1_public.pop("_h6", None)

    (out / "GB1_INVARIANTS.json").write_text(
        json.dumps(gb1_public["invariants"], indent=2), encoding="utf-8"
    )
    (out / "GB1_DIAGNOSTICS.json").write_text(
        json.dumps(gb1_public, indent=2), encoding="utf-8"
    )
    (out / "H6_RESIDUAL_STRUCTURE.json").write_text(
        json.dumps(h6, indent=2), encoding="utf-8"
    )
    (out / "STAGE2_MATRIX.json").write_text(
        json.dumps(matrix, indent=2), encoding="utf-8"
    )

    manifest = {
        "version": "NABU_ROOT_CAUSE_FALSIFICATION_STAGE2_MANIFEST_V1",
        "trpb_sha256": sha256_file(Path(args.trpb)),
        "gb1_sha256": sha256_file(Path(args.gb1)),
        "gb1_git_blob_sha": git_blob_sha(Path(args.gb1)),
        "gb1_source_commit": "0861a3f1f6c58a9597fde689eb7f9a07b38e2638",
        "gb1_wt": GB1_WT,
        "seed": SEED,
        "nucb_consumed": False,
        "phase3_opened": False,
    }
    (out / "RUN_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )

    hashes = {}
    for name in (
        "GB1_INVARIANTS.json",
        "GB1_DIAGNOSTICS.json",
        "H6_RESIDUAL_STRUCTURE.json",
        "STAGE2_MATRIX.json",
        "RUN_MANIFEST.json",
    ):
        hashes[name] = sha256_file(out / name)
    (out / "OUTPUT_HASHES.json").write_text(
        json.dumps(hashes, indent=2), encoding="utf-8"
    )

    print(json.dumps(matrix, indent=2))


if __name__ == "__main__":
    main()
