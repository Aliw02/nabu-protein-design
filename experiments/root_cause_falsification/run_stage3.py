from __future__ import annotations

import argparse
import hashlib
import json
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import run_falsification as s1
import run_stage2 as s2

from nabu_protein.v83 import NabuV83Model
from nabu_protein.higher_order import additive_score
from run_v9_pair_transfer_ired import derive_reference, mutation_set

SEED = 161
EPS = 1e-12


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_trpb(path: Path):
    fit, _, test, source_info = s1.split_frame(path, "TrpB")
    fit_seq = fit["sequence"].astype(str).tolist()
    test_seq = test["sequence"].astype(str).tolist()
    ref = derive_reference(fit_seq)
    fit_sets = [mutation_set(x, ref) for x in fit_seq]
    test_sets = [mutation_set(x, ref) for x in test_seq]
    fit_y = fit["target"].to_numpy(dtype=float)
    test_y = test["target"].to_numpy(dtype=float)
    orders = np.asarray([len(x) for x in test_sets], dtype=int)
    ids = test["row_id"].astype(str).tolist()
    return {
        "name": "TrpB",
        "fit_ids": fit["row_id"].astype(str).tolist(),
        "test_ids": ids,
        "fit_sets": fit_sets,
        "test_sets": test_sets,
        "fit_y": fit_y,
        "test_y": test_y,
        "orders": orders,
        "source_info": source_info,
    }


def load_gb1(path: Path):
    if path.stat().st_size != s2.GB1_EXPECTED_SIZE:
        raise RuntimeError("GB1 source size drift.")
    if s2.git_blob_sha(path) != s2.GB1_GIT_BLOB_SHA:
        raise RuntimeError("GB1 git blob identity drift.")

    frame = pd.read_csv(path)
    frame = frame[["variant", "fitness"]].dropna().copy()
    frame["variant"] = frame["variant"].astype(str).str.strip().str.upper()
    frame["fitness"] = pd.to_numeric(frame["fitness"], errors="raise")
    frame = frame[
        frame["variant"].str.len().eq(4)
        & frame["variant"].map(lambda x: set(x).issubset(s2.STANDARD_AA))
    ].copy()
    if frame["variant"].duplicated().any():
        raise RuntimeError("GB1 duplicate variants.")
    if s2.GB1_WT not in set(frame["variant"]):
        raise RuntimeError("GB1 WT missing.")

    frame["mutation_set"] = frame["variant"].map(s2.gb1_mutation_set)
    frame["mutation_order"] = frame["mutation_set"].map(len)
    fit = frame[frame["mutation_order"].le(2)].copy()
    test = frame[frame["mutation_order"].ge(3)].copy()

    return {
        "name": "GB1",
        "fit_ids": fit["variant"].astype(str).tolist(),
        "test_ids": test["variant"].astype(str).tolist(),
        "fit_sets": fit["mutation_set"].tolist(),
        "test_sets": test["mutation_set"].tolist(),
        "fit_y": fit["fitness"].to_numpy(dtype=float),
        "test_y": test["fitness"].to_numpy(dtype=float),
        "orders": test["mutation_order"].to_numpy(dtype=int),
        "source_info": {
            "git_blob_sha": s2.git_blob_sha(path),
            "sha256": sha256_file(path),
        },
    }


class IdentityTransform:
    name = "identity"

    def fit(self, y):
        return self

    def transform(self, y):
        return np.asarray(y, dtype=float)


def fit_transform(name: str, fit_y):
    if name == "identity":
        t = IdentityTransform()
    elif name == "robust_asinh":
        t = s1.RobustAsinhTransform()
    elif name == "yeo_johnson":
        t = s1.YeoJohnsonTransform()
    else:
        raise ValueError(name)
    t.fit(fit_y)
    return t


def exact_scale_decomposition(dataset, transform_name: str):
    transform = fit_transform(transform_name, dataset["fit_y"])
    fit_y = transform.transform(dataset["fit_y"])
    test_y = transform.transform(dataset["test_y"])

    lookup, duplicates = s1.unique_lookup(
        dataset["fit_sets"], fit_y, f"{dataset['name']}_{transform_name}_fit"
    )
    if duplicates:
        raise RuntimeError(
            f"{dataset['name']} {transform_name}: duplicate low-order mutation sets."
        )
    wt, e1, e2 = s1.exact_components(lookup)

    keys = sorted(e2)
    values = np.asarray([e2[k] for k in keys], dtype=float)
    shuffled_values = np.random.default_rng(SEED).permutation(values)
    shuffled_e2 = {k: float(v) for k, v in zip(keys, shuffled_values)}

    o1, o2, shuffled, _, o2_support, _ = s1.o1_o2_predictions(
        dataset["test_sets"], wt, e1, e2, shuffled_e2
    )

    by_order = {}
    for order in sorted(set(dataset["orders"].tolist())):
        mask = (dataset["orders"] == order) & o2_support
        if int(np.sum(mask)) >= 2:
            by_order[str(order)] = {
                "support_count": int(np.sum(mask)),
                "O1": s1.masked_metrics(test_y, o1, mask),
                "O2": s1.masked_metrics(test_y, o2, mask),
                "SHUFFLED_O2": s1.masked_metrics(test_y, shuffled, mask),
            }

    # Triple e3 and context diagnostics.
    triple_e3 = {}
    base_pairs = []
    conditional_pairs = []
    triple_o2_values = []
    triple_e3_values = []
    triple_targets = []

    for muts, y, order, supported, pred_o2 in zip(
        dataset["test_sets"],
        test_y,
        dataset["orders"],
        o2_support,
        o2,
    ):
        if order != 3 or not supported:
            continue
        e3 = float(y - pred_o2)
        triple_e3[tuple(muts)] = e3
        triple_o2_values.append(float(pred_o2))
        triple_e3_values.append(e3)
        triple_targets.append(float(y))

        for pair in combinations(muts, 2):
            pair = tuple(pair)
            base = float(e2[pair])
            base_pairs.append(base)
            conditional_pairs.append(float(base + e3))

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

    triple_e3_values = np.asarray(triple_e3_values, dtype=float)
    triple_o2_values = np.asarray(triple_o2_values, dtype=float)
    triple_targets = np.asarray(triple_targets, dtype=float)

    target_sd = float(np.std(triple_targets)) if len(triple_targets) else 0.0
    e3_rms = (
        None
        if len(triple_e3_values) == 0
        else float(np.sqrt(np.mean(triple_e3_values ** 2)))
    )

    context = {
        "supported_triple_count": int(len(triple_e3)),
        "pair_context_instance_count": int(len(base_pairs)),
        "sign_flip_rate": sign_flip,
        "spearman_base_vs_conditional": s1.finite_spearman(
            base_pairs, conditional_pairs
        ),
        "pearson_base_vs_conditional": s1.finite_pearson(
            base_pairs, conditional_pairs
        ),
        "e3_rms": e3_rms,
        "triple_target_std": target_sd,
        "e3_rms_over_target_std": (
            None
            if e3_rms is None or target_sd <= EPS
            else float(e3_rms / target_sd)
        ),
        "spearman_e3_vs_o2": s1.finite_spearman(
            triple_e3_values, triple_o2_values
        ),
        "pearson_e3_vs_o2": s1.finite_pearson(
            triple_e3_values, triple_o2_values
        ),
    }

    # Quad decomposition using measured constituent triple e3 values.
    o3 = np.full(len(test_y), np.nan, dtype=float)
    e4 = np.full(len(test_y), np.nan, dtype=float)
    o3_support = np.zeros(len(test_y), dtype=bool)

    for i, (muts, y, order, supported) in enumerate(
        zip(dataset["test_sets"], test_y, dataset["orders"], o2_support)
    ):
        if order != 4 or not supported:
            continue
        triples = [tuple(x) for x in combinations(muts, 3)]
        if all(t in triple_e3 for t in triples):
            o3_support[i] = True
            o3[i] = float(o2[i] + sum(triple_e3[t] for t in triples))
            e4[i] = float(y - o3[i])

    quad_mask = o3_support & (dataset["orders"] == 4)
    if int(np.sum(quad_mask)) >= 2:
        e4_vals = e4[quad_mask]
        quad_target = test_y[quad_mask]
        quad_sd = float(np.std(quad_target))
        e4_rms = float(np.sqrt(np.mean(e4_vals ** 2)))
        quad = {
            "support_count": int(np.sum(quad_mask)),
            "O2": s1.masked_metrics(test_y, o2, quad_mask),
            "O3": s1.masked_metrics(test_y, o3, quad_mask),
            "O3_minus_O2_spearman": float(
                s1.masked_metrics(test_y, o3, quad_mask)["spearman"]
                - s1.masked_metrics(test_y, o2, quad_mask)["spearman"]
            ),
            "e4_rms": e4_rms,
            "quad_target_std": quad_sd,
            "e4_rms_over_target_std": (
                None if quad_sd <= EPS else float(e4_rms / quad_sd)
            ),
            "spearman_e4_vs_o3": s1.finite_spearman(
                e4_vals, o3[quad_mask]
            ),
            "pearson_e4_vs_o3": s1.finite_pearson(
                e4_vals, o3[quad_mask]
            ),
        }
    else:
        quad = {"support_count": int(np.sum(quad_mask))}

    meta = {}
    if transform_name == "robust_asinh":
        meta = {
            "center": float(transform.center),
            "scale": float(transform.scale),
            "scale_source": transform.scale_source,
        }
    elif transform_name == "yeo_johnson":
        meta = {"lambda": float(transform.pt.lambdas_[0])}

    return {
        "transform": transform_name,
        "fit_only_parameters": meta,
        "by_order": by_order,
        "context": context,
        "quad_decomposition": quad,
        "_internal": {
            "test_target": test_y,
            "o2": o2,
            "o2_support": o2_support,
            "triple_e3": triple_e3,
            "wt": wt,
            "e1": e1,
            "e2": e2,
        },
    }


def train_b2(dataset):
    model = NabuV83Model().fit(
        mutation_sets=dataset["fit_sets"],
        labels=dataset["fit_y"],
        candidate_ids=dataset["fit_ids"],
    )
    base = model.model["base"]
    mean = float(base["global_mean"])
    b2 = np.empty(len(dataset["test_sets"]), dtype=float)
    for i, muts in enumerate(dataset["test_sets"]):
        if all(m in base["main"] for m in muts):
            b2[i] = float(additive_score(muts, mean, base["main"]))
        else:
            b2[i] = mean
    return b2


def affine_fit(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    A = np.column_stack([np.ones(len(x), dtype=float), x])
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    return float(beta[0]), float(beta[1])


def affine_predict(x, a, b):
    return a + b * np.asarray(x, dtype=float)


def cross_order_transfer(dataset, identity_diag):
    internal = identity_diag["_internal"]
    y = internal["test_target"]
    o2 = internal["o2"]
    support = internal["o2_support"]
    orders = dataset["orders"]

    triple_mask = (orders == 3) & support
    quad_mask = (orders == 4) & support
    if int(np.sum(triple_mask)) < 100 or int(np.sum(quad_mask)) < 100:
        return {
            "evaluable": False,
            "triple_count": int(np.sum(triple_mask)),
            "quad_count": int(np.sum(quad_mask)),
        }

    residual3 = y[triple_mask] - o2[triple_mask]
    x3 = o2[triple_mask]
    a, b = affine_fit(x3, residual3)

    rng = np.random.default_rng(SEED)
    shuffled_residual3 = rng.permutation(residual3)
    sa, sb = affine_fit(x3, shuffled_residual3)

    quad_o2 = o2[quad_mask]
    pred = quad_o2 + affine_predict(quad_o2, a, b)
    shuffled_pred = quad_o2 + affine_predict(quad_o2, sa, sb)
    b2 = train_b2(dataset)[quad_mask]
    target = y[quad_mask]

    return {
        "evaluable": True,
        "triple_count": int(np.sum(triple_mask)),
        "quad_count": int(np.sum(quad_mask)),
        "triple_affine": {"intercept": a, "slope": b},
        "shuffled_triple_affine": {"intercept": sa, "slope": sb},
        "quad_metrics": {
            "RAW_O2": s1.extended_metrics(target, quad_o2),
            "TRIPLE_AFFINE_LAW": s1.extended_metrics(target, pred),
            "SHUFFLED_TRIPLE_LAW": s1.extended_metrics(
                target, shuffled_pred
            ),
            "B2": s1.extended_metrics(target, b2),
        },
    }


def public_diag(diag):
    out = dict(diag)
    out.pop("_internal", None)
    return out


def classify(trpb_scales, gb1_scales, transfer):
    identity_t = trpb_scales["identity"]
    robust_t = trpb_scales["robust_asinh"]
    identity_g = gb1_scales["identity"]
    robust_g = gb1_scales["robust_asinh"]

    # Diagnostic A: same-order scale rescue.
    scale_rescue_orders = []
    common_orders = sorted(
        set(identity_t["by_order"])
        & set(robust_t["by_order"])
        & set(identity_g["by_order"])
        & set(robust_g["by_order"]),
        key=int,
    )
    for order in common_orders:
        t_gain = float(
            robust_t["by_order"][order]["O2"]["spearman"]
            - identity_t["by_order"][order]["O2"]["spearman"]
        )
        g_gain = float(
            robust_g["by_order"][order]["O2"]["spearman"]
            - identity_g["by_order"][order]["O2"]["spearman"]
        )
        scale_rescue_orders.append(
            {
                "order": order,
                "TrpB_O2_spearman_gain": t_gain,
                "GB1_O2_spearman_gain": g_gain,
                "primary_scale_rescue": bool(
                    t_gain >= 0.10 and g_gain >= 0.10
                ),
            }
        )
    scale_rescue = any(x["primary_scale_rescue"] for x in scale_rescue_orders)

    # Diagnostic B.
    def reductions(identity, robust):
        i = identity["context"]
        r = robust["context"]
        sf_red = (
            None
            if i["sign_flip_rate"] is None
            or r["sign_flip_rate"] is None
            or i["sign_flip_rate"] <= EPS
            else float(1.0 - r["sign_flip_rate"] / i["sign_flip_rate"])
        )
        e3_red = (
            None
            if i["e3_rms_over_target_std"] is None
            or r["e3_rms_over_target_std"] is None
            or i["e3_rms_over_target_std"] <= EPS
            else float(
                1.0
                - r["e3_rms_over_target_std"]
                / i["e3_rms_over_target_std"]
            )
        )
        mediated = bool(
            sf_red is not None
            and e3_red is not None
            and sf_red >= 0.50
            and e3_red >= 0.50
        )
        persistent = bool(
            r["sign_flip_rate"] is not None
            and r["e3_rms_over_target_std"] is not None
            and r["sign_flip_rate"] >= 0.20
            and r["e3_rms_over_target_std"] >= 0.50
        )
        return {
            "sign_flip_reduction": sf_red,
            "e3_ratio_reduction": e3_red,
            "scale_mediated": mediated,
            "persistent_after_scale": persistent,
        }

    red_t = reductions(identity_t, robust_t)
    red_g = reductions(identity_g, robust_g)
    mediated_both = red_t["scale_mediated"] and red_g["scale_mediated"]
    persistent_both = (
        red_t["persistent_after_scale"]
        and red_g["persistent_after_scale"]
    )
    if mediated_both and not persistent_both:
        context_outcome = "SCALE_MEDIATES_MOST_CONTEXT_DEPENDENCE"
    elif persistent_both and not mediated_both:
        context_outcome = "CONTEXT_DEPENDENCE_PERSISTS_AFTER_SCALE_CORRECTION"
    else:
        context_outcome = "MIXED_OR_UNRESOLVED"

    # Diagnostic C.
    transfer_by_dataset = {}
    transfer_supported = True
    for name in ("TrpB", "GB1"):
        block = transfer[name]
        if not block.get("evaluable"):
            transfer_by_dataset[name] = {
                "supported": False,
                "reason": "not_evaluable",
            }
            transfer_supported = False
            continue
        q = block["quad_metrics"]
        real = q["TRIPLE_AFFINE_LAW"]["spearman"]
        raw = q["RAW_O2"]["spearman"]
        shuffled = q["SHUFFLED_TRIPLE_LAW"]["spearman"]
        ok = bool(
            real >= raw + 0.10
            and real >= shuffled + 0.10
        )
        transfer_by_dataset[name] = {
            "supported": ok,
            "triple_law_minus_raw_o2": float(real - raw),
            "triple_law_minus_shuffled": float(real - shuffled),
        }
        transfer_supported = transfer_supported and ok

    # Diagnostic D.
    robust_e4 = {
        "TrpB": robust_t["quad_decomposition"].get(
            "e4_rms_over_target_std"
        ),
        "GB1": robust_g["quad_decomposition"].get(
            "e4_rms_over_target_std"
        ),
    }
    robust_e4_nonmaterial = all(
        value is not None and value < 0.50
        for value in robust_e4.values()
    )

    promote = bool(
        scale_rescue
        and context_outcome == "SCALE_MEDIATES_MOST_CONTEXT_DEPENDENCE"
        and transfer_supported
        and robust_e4_nonmaterial
    )

    root_status = (
        "OBSERVATION_SCALE_COMPOSITION_DOMINANT_MECHANISM"
        if promote
        else "ROOT_CAUSE_NOT_YET_IDENTIFIED"
    )

    return {
        "version": "NABU_ROOT_CAUSE_FALSIFICATION_STAGE3_V1",
        "root_cause_status": root_status,
        "diagnostic_A_scale_rescue": {
            "passes_primary_gate": scale_rescue,
            "by_order": scale_rescue_orders,
        },
        "diagnostic_B_scale_vs_context": {
            "outcome": context_outcome,
            "TrpB": red_t,
            "GB1": red_g,
        },
        "diagnostic_C_cross_order_transfer": {
            "supported_on_both": transfer_supported,
            "by_dataset": transfer_by_dataset,
        },
        "diagnostic_D_robust_e4": {
            "e4_rms_over_target_std": robust_e4,
            "nonmaterial_on_both": robust_e4_nonmaterial,
        },
        "promotion_gate": {
            "scale_rescue": scale_rescue,
            "scale_mediates_context": (
                context_outcome
                == "SCALE_MEDIATES_MOST_CONTEXT_DEPENDENCE"
            ),
            "cross_order_transfer": transfer_supported,
            "robust_e4_nonmaterial": robust_e4_nonmaterial,
            "all_required": promote,
        },
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

    datasets = {
        "TrpB": load_trpb(Path(args.trpb)),
        "GB1": load_gb1(Path(args.gb1)),
    }

    scales = {}
    for name, ds in datasets.items():
        scales[name] = {
            "identity": exact_scale_decomposition(ds, "identity"),
            "robust_asinh": exact_scale_decomposition(
                ds, "robust_asinh"
            ),
            "yeo_johnson": exact_scale_decomposition(
                ds, "yeo_johnson"
            ),
        }

    transfer = {
        name: cross_order_transfer(ds, scales[name]["identity"])
        for name, ds in datasets.items()
    }

    matrix = classify(
        scales["TrpB"],
        scales["GB1"],
        transfer,
    )

    scale_public = {
        dataset: {
            transform: public_diag(diag)
            for transform, diag in blocks.items()
        }
        for dataset, blocks in scales.items()
    }
    quad = {
        dataset: {
            transform: blocks[transform]["quad_decomposition"]
            for transform in ("identity", "robust_asinh", "yeo_johnson")
        }
        for dataset, blocks in scales.items()
    }

    (out / "STAGE3_MATRIX.json").write_text(
        json.dumps(matrix, indent=2), encoding="utf-8"
    )
    (out / "SCALE_CONTEXT_DIAGNOSTICS.json").write_text(
        json.dumps(scale_public, indent=2), encoding="utf-8"
    )
    (out / "CROSS_ORDER_TRANSFER.json").write_text(
        json.dumps(transfer, indent=2), encoding="utf-8"
    )
    (out / "QUAD_DECOMPOSITION.json").write_text(
        json.dumps(quad, indent=2), encoding="utf-8"
    )

    manifest = {
        "version": "NABU_ROOT_CAUSE_STAGE3_RUN_MANIFEST_V1",
        "seed": SEED,
        "trpb_sha256": sha256_file(Path(args.trpb)),
        "gb1_sha256": sha256_file(Path(args.gb1)),
        "gb1_git_blob_sha": s2.git_blob_sha(Path(args.gb1)),
        "nucb_consumed": False,
        "phase3_opened": False,
    }
    (out / "RUN_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )

    hashes = {}
    for name in (
        "STAGE3_MATRIX.json",
        "SCALE_CONTEXT_DIAGNOSTICS.json",
        "CROSS_ORDER_TRANSFER.json",
        "QUAD_DECOMPOSITION.json",
        "RUN_MANIFEST.json",
    ):
        hashes[name] = sha256_file(out / name)
    (out / "OUTPUT_HASHES.json").write_text(
        json.dumps(hashes, indent=2), encoding="utf-8"
    )

    print(json.dumps(matrix, indent=2))


if __name__ == "__main__":
    main()
