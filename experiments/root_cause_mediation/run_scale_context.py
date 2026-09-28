from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
FALSIFICATION = HERE.parent / "root_cause_falsification"
sys.path.insert(0, str(FALSIFICATION))
DEV = HERE.parent / "v9_pair_transfer_dev"
sys.path.insert(0, str(DEV))

from run_falsification import (
    EPS,
    SEED,
    IdentityTransform,
    RobustAsinhTransform,
    exact_components,
    masked_metrics,
    o1_o2_predictions,
    pair_context_diagnostics,
    split_frame,
    unique_lookup,
)
from run_stage2 import (
    GB1_EXPECTED_SIZE,
    GB1_GIT_BLOB_SHA,
    GB1_WT,
    gb1_mutation_set,
    git_blob_sha,
)
from run_v9_pair_transfer_ired import derive_reference, mutation_set


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_trpb(path: Path):
    fit, _, test, source = split_frame(path, "TrpB")
    reference = derive_reference(fit["sequence"].astype(str).tolist())
    fit_sets = [
        mutation_set(seq, reference)
        for seq in fit["sequence"].astype(str)
    ]
    test_sets = [
        mutation_set(seq, reference)
        for seq in test["sequence"].astype(str)
    ]
    return {
        "name": "TrpB",
        "fit_sets": fit_sets,
        "fit_target": fit["target"].to_numpy(dtype=float),
        "test_sets": test_sets,
        "test_target": test["target"].to_numpy(dtype=float),
        "test_orders": np.asarray([len(x) for x in test_sets], dtype=int),
        "source": source,
        "reference": reference,
    }


def load_gb1(path: Path):
    if path.stat().st_size != GB1_EXPECTED_SIZE:
        raise RuntimeError("GB1 source size drift.")
    if git_blob_sha(path) != GB1_GIT_BLOB_SHA:
        raise RuntimeError("GB1 source Git blob drift.")

    frame = pd.read_csv(path)
    required = {"variant", "fitness"}
    missing = required - set(frame.columns)
    if missing:
        raise RuntimeError(f"GB1 missing columns: {sorted(missing)}")

    frame = frame[["variant", "fitness"]].dropna().copy()
    frame["variant"] = frame["variant"].astype(str).str.strip().str.upper()
    frame["fitness"] = pd.to_numeric(frame["fitness"], errors="raise")
    frame["mutation_set"] = frame["variant"].map(gb1_mutation_set)
    frame["mutation_order"] = frame["mutation_set"].map(len)

    if frame["variant"].duplicated().any():
        raise RuntimeError("GB1 duplicate variant.")
    if GB1_WT not in set(frame["variant"]):
        raise RuntimeError("GB1 WT missing.")

    wt_value = float(
        frame.loc[frame["variant"].eq(GB1_WT), "fitness"].iloc[0]
    )
    if abs(wt_value - 1.0) > 1e-12:
        raise RuntimeError(f"GB1 WT fitness drift: {wt_value}")

    fit = frame[frame["mutation_order"].le(2)].copy()
    test = frame[frame["mutation_order"].ge(3)].copy()

    return {
        "name": "GB1",
        "fit_sets": fit["mutation_set"].tolist(),
        "fit_target": fit["fitness"].to_numpy(dtype=float),
        "test_sets": test["mutation_set"].tolist(),
        "test_target": test["fitness"].to_numpy(dtype=float),
        "test_orders": test["mutation_order"].to_numpy(dtype=int),
        "source": {
            "git_blob_sha": GB1_GIT_BLOB_SHA,
            "size": GB1_EXPECTED_SIZE,
            "wt": GB1_WT,
        },
        "reference": GB1_WT,
    }


def portability_block(
    target,
    orders,
    o1,
    o2,
    shuffled,
    support,
):
    result = {}
    for order in sorted(set(orders.tolist())):
        mask = (orders == order) & support
        count = int(np.sum(mask))
        if count < 100:
            result[str(order)] = {
                "evaluable": False,
                "support_count": count,
            }
            continue

        m1 = masked_metrics(target, o1, mask)
        m2 = masked_metrics(target, o2, mask)
        msh = masked_metrics(target, shuffled, mask)
        s1 = m1["spearman"]
        s2 = m2["spearman"]
        ssh = msh["spearman"]
        if s1 is None or s2 is None or ssh is None:
            result[str(order)] = {
                "evaluable": False,
                "support_count": count,
            }
            continue

        result[str(order)] = {
            "evaluable": True,
            "support_count": count,
            "O1": m1,
            "O2": m2,
            "SHUFFLED_O2": msh,
            "o2_minus_o1": float(s2 - s1),
            "o2_minus_shuffled": float(s2 - ssh),
            "non_portable": bool(
                (s2 - s1 < 0.05) or (s2 <= ssh + 0.02)
            ),
            "portable": bool(
                (s2 - s1 >= 0.05) and (s2 > ssh + 0.02)
            ),
        }
    return result


def analyze_scale(data, transform):
    transform.fit(data["fit_target"])
    fit_y = transform.transform(data["fit_target"])
    test_y = transform.transform(data["test_target"])

    lookup, duplicates = unique_lookup(
        data["fit_sets"], fit_y, f"{data['name']}_{transform.name}"
    )
    if duplicates:
        raise RuntimeError(
            f"{data['name']} duplicate low-order mutation sets."
        )

    wt, e1, e2 = exact_components(lookup)
    pair_keys = sorted(e2)
    pair_values = np.asarray([e2[k] for k in pair_keys], dtype=float)
    shuffled_values = np.random.default_rng(SEED).permutation(pair_values)
    shuffled_e2 = {
        key: float(value)
        for key, value in zip(pair_keys, shuffled_values)
    }

    o1, o2, shuffled, _, support, _ = o1_o2_predictions(
        data["test_sets"],
        wt,
        e1,
        e2,
        shuffled_e2,
    )
    portability = portability_block(
        test_y,
        data["test_orders"],
        o1,
        o2,
        shuffled,
        support,
    )
    context = pair_context_diagnostics(
        data["test_sets"],
        test_y,
        data["test_orders"],
        wt,
        e1,
        e2,
    )

    h3_order3 = portability.get("3", {})
    h4_material = bool(
        context["supported_triple_count"] >= 100
        and context["sign_flip_rate"] is not None
        and context["e3_rms_over_target_std"] is not None
        and context["sign_flip_rate"] >= 0.20
        and context["e3_rms_over_target_std"] >= 0.50
    )

    parameters = {}
    if isinstance(transform, RobustAsinhTransform):
        parameters = {
            "center": float(transform.center),
            "scale": float(transform.scale),
            "scale_source": transform.scale_source,
        }

    return {
        "scale": transform.name,
        "fit_only_parameters": parameters,
        "exact_single_count": int(len(e1)),
        "exact_pair_count": int(len(e2)),
        "portability_by_order": portability,
        "conditional_pair_context": context,
        "H3_order3_non_portable": bool(
            h3_order3.get("evaluable")
            and h3_order3.get("non_portable")
        ),
        "H3_order3_portable": bool(
            h3_order3.get("evaluable")
            and h3_order3.get("portable")
        ),
        "H4_material": h4_material,
    }


def dataset_analysis(data):
    identity = analyze_scale(data, IdentityTransform())
    robust = analyze_scale(data, RobustAsinhTransform())

    def ratio(new, old):
        if new is None or old is None or abs(old) <= EPS:
            return None
        return float(new / old)

    identity_context = identity["conditional_pair_context"]
    robust_context = robust["conditional_pair_context"]

    mediation = {
        "sign_flip_ratio_robust_over_identity": ratio(
            robust_context["sign_flip_rate"],
            identity_context["sign_flip_rate"],
        ),
        "e3_rms_over_target_std_ratio_robust_over_identity": ratio(
            robust_context["e3_rms_over_target_std"],
            identity_context["e3_rms_over_target_std"],
        ),
        "conditional_pair_spearman_change": (
            None
            if robust_context["spearman_base_vs_conditional"] is None
            or identity_context["spearman_base_vs_conditional"] is None
            else float(
                robust_context["spearman_base_vs_conditional"]
                - identity_context["spearman_base_vs_conditional"]
            )
        ),
    }

    return {
        "dataset": data["name"],
        "source": data["source"],
        "identity": identity,
        "robust_asinh": robust,
        "mediation_descriptives": mediation,
    }


def decision(trpb, gb1):
    robusts = [trpb["robust_asinh"], gb1["robust_asinh"]]

    sufficient = all(
        item["H3_order3_portable"] and not item["H4_material"]
        for item in robusts
    )
    persists = all(
        item["H3_order3_non_portable"] and item["H4_material"]
        for item in robusts
    )

    if sufficient:
        label = "SCALE_CORRECTION_SUFFICIENT_FOR_H3_H4"
    elif persists:
        label = "CONTEXT_DEPENDENCE_PERSISTS_AFTER_STANDARD_SCALE_CORRECTION"
    else:
        label = "SCALE_VS_CONTEXT_UNRESOLVED"

    return {
        "version": "NABU_SCALE_CONTEXT_MEDIATION_V1",
        "decision": label,
        "root_cause_status": "ROOT_CAUSE_NOT_YET_IDENTIFIED",
        "TrpB_robust_H3_non_portable": trpb["robust_asinh"][
            "H3_order3_non_portable"
        ],
        "TrpB_robust_H4_material": trpb["robust_asinh"]["H4_material"],
        "GB1_robust_H3_non_portable": gb1["robust_asinh"][
            "H3_order3_non_portable"
        ],
        "GB1_robust_H4_material": gb1["robust_asinh"]["H4_material"],
        "interpretation": (
            "This tests whether the prespecified robust_asinh low-order-only "
            "scale correction is sufficient to remove pair non-portability and "
            "conditional pair context dependence. It does not adjudicate every "
            "possible nonlinear latent mapping."
        ),
        "nucb_consumed": False,
        "phase3_opened": False,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--trpb", required=True)
    parser.add_argument("--gb1", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    trpb_path = Path(args.trpb)
    gb1_path = Path(args.gb1)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    trpb = dataset_analysis(load_trpb(trpb_path))
    gb1 = dataset_analysis(load_gb1(gb1_path))
    matrix = decision(trpb, gb1)

    (out / "TRPB_SCALE_CONTEXT.json").write_text(
        json.dumps(trpb, indent=2), encoding="utf-8"
    )
    (out / "GB1_SCALE_CONTEXT.json").write_text(
        json.dumps(gb1, indent=2), encoding="utf-8"
    )
    (out / "SCALE_CONTEXT_MATRIX.json").write_text(
        json.dumps(matrix, indent=2), encoding="utf-8"
    )

    manifest = {
        "version": "NABU_SCALE_CONTEXT_MEDIATION_RUN_V1",
        "seed": SEED,
        "trpb_sha256": sha256_file(trpb_path),
        "gb1_sha256": sha256_file(gb1_path),
        "gb1_git_blob_sha": git_blob_sha(gb1_path),
        "nucb_consumed": False,
        "phase3_opened": False,
    }
    (out / "RUN_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )

    hashes = {}
    for name in (
        "TRPB_SCALE_CONTEXT.json",
        "GB1_SCALE_CONTEXT.json",
        "SCALE_CONTEXT_MATRIX.json",
        "RUN_MANIFEST.json",
    ):
        hashes[name] = sha256_file(out / name)
    (out / "OUTPUT_HASHES.json").write_text(
        json.dumps(hashes, indent=2), encoding="utf-8"
    )

    print(json.dumps(matrix, indent=2))


if __name__ == "__main__":
    main()
