from __future__ import annotations

import argparse
import hashlib
import json
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.stats import pearsonr, spearmanr
from sklearn.ensemble import RandomForestRegressor

HERE = Path(__file__).resolve().parent
FALS = HERE.parent / "root_cause_falsification"
sys.path.insert(0, str(FALS))

import run_falsification as s1
import run_stage3 as s3

SEED = 161
EPS = 1e-12


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
    x = float(spearmanr(a, b).statistic)
    return x if np.isfinite(x) else None


def finite_pearson(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) < 2 or np.std(a) <= EPS or np.std(b) <= EPS:
        return None
    x = float(pearsonr(a, b).statistic)
    return x if np.isfinite(x) else None


def build_triples(dataset):
    transform = s1.RobustAsinhTransform()
    transform.fit(dataset["fit_y"])
    fit_y = transform.transform(dataset["fit_y"])
    test_y = transform.transform(dataset["test_y"])

    lookup, duplicates = s1.unique_lookup(
        dataset["fit_sets"], fit_y, f"{dataset['name']}_prior_fit"
    )
    if duplicates:
        raise RuntimeError(f"{dataset['name']}: duplicate low-order mutation sets.")

    wt, e1, e2 = s1.exact_components(lookup)
    sigma_low = float(np.std(fit_y))
    if sigma_low <= EPS:
        raise RuntimeError(f"{dataset['name']}: zero low-order scale.")

    X = []
    y = []
    o2_values = []
    ids = []

    for cid, muts, target, order in zip(
        dataset["test_ids"],
        dataset["test_sets"],
        test_y,
        dataset["orders"],
    ):
        if int(order) != 3:
            continue
        if not all(m in e1 for m in muts):
            continue
        pairs = [tuple(pair) for pair in combinations(muts, 2)]
        if not all(pair in e2 for pair in pairs):
            continue

        e1v = np.asarray([e1[m] for m in muts], dtype=float)
        e2v = np.asarray([e2[p] for p in pairs], dtype=float)
        o1 = float(wt + np.sum(e1v))
        o2 = float(o1 + np.sum(e2v))
        e3 = float(target - o2)

        X.append(
            [
                o1 / sigma_low,
                o2 / sigma_low,
                float(np.sum(e1v)) / sigma_low,
                float(np.mean(np.abs(e1v))) / sigma_low,
                float(np.max(np.abs(e1v))) / sigma_low,
                float(np.std(e1v)) / sigma_low,
                float(np.min(e1v)) / sigma_low,
                float(np.max(e1v)) / sigma_low,
                float(np.sum(e2v)) / sigma_low,
                float(np.mean(np.abs(e2v))) / sigma_low,
                float(np.max(np.abs(e2v))) / sigma_low,
                float(np.std(e2v)) / sigma_low,
                float(np.min(e2v)) / sigma_low,
                float(np.max(e2v)) / sigma_low,
                float(np.mean(e2v > 0.0)),
            ]
        )
        y.append(e3 / sigma_low)
        o2_values.append(o2 / sigma_low)
        ids.append(str(cid))

    return {
        "dataset": dataset["name"],
        "X": np.asarray(X, dtype=float),
        "y": np.asarray(y, dtype=float),
        "o2": np.asarray(o2_values, dtype=float),
        "ids": np.asarray(ids, dtype=str),
        "sigma_low": sigma_low,
        "transform": {
            "center": float(transform.center),
            "scale": float(transform.scale),
            "scale_source": transform.scale_source,
        },
    }


def residualize_against_o2(o2, y):
    o2 = np.asarray(o2, dtype=float)
    y = np.asarray(y, dtype=float)
    A = np.column_stack([np.ones(len(o2), dtype=float), o2])
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    fitted = A @ beta
    return y - fitted, {
        "intercept": float(beta[0]),
        "slope": float(beta[1]),
    }


def transfer(source, target):
    if len(source["y"]) < 100 or len(target["y"]) < 100:
        raise RuntimeError("Insufficient supported triples for prior transfer.")

    model = RandomForestRegressor(
        n_estimators=300,
        min_samples_leaf=50,
        max_features="sqrt",
        random_state=SEED,
        n_jobs=1,
    )
    model.fit(source["X"], source["y"])
    pred = np.asarray(model.predict(target["X"]), dtype=float)

    rng = np.random.default_rng(SEED)
    shuffled_source_y = rng.permutation(source["y"])
    control = RandomForestRegressor(
        n_estimators=300,
        min_samples_leaf=50,
        max_features="sqrt",
        random_state=SEED,
        n_jobs=1,
    )
    control.fit(source["X"], shuffled_source_y)
    shuffled_pred = np.asarray(control.predict(target["X"]), dtype=float)

    real_s = finite_spearman(target["y"], pred)
    shuffled_s = finite_spearman(target["y"], shuffled_pred)

    true_conditional, true_affine = residualize_against_o2(
        target["o2"], target["y"]
    )
    pred_conditional, pred_affine = residualize_against_o2(
        target["o2"], pred
    )

    return {
        "source": source["dataset"],
        "target": target["dataset"],
        "source_count": int(len(source["y"])),
        "target_count": int(len(target["y"])),
        "target_spearman": real_s,
        "target_pearson": finite_pearson(target["y"], pred),
        "target_rmse": float(
            np.sqrt(np.mean((target["y"] - pred) ** 2))
        ),
        "shuffled_control_spearman": shuffled_s,
        "real_minus_shuffled_spearman": (
            None
            if real_s is None or shuffled_s is None
            else float(real_s - shuffled_s)
        ),
        "conditional_residual_spearman": finite_spearman(
            true_conditional, pred_conditional
        ),
        "target_true_e3_vs_o2_affine": true_affine,
        "prediction_vs_o2_affine": pred_affine,
        "source_sigma_low": float(source["sigma_low"]),
        "target_sigma_low": float(target["sigma_low"]),
        "source_transform": source["transform"],
        "target_transform": target["transform"],
    }


def classify(gb_to_tr, tr_to_gb):
    blocks = [gb_to_tr, tr_to_gb]

    supported = all(
        b["target_spearman"] is not None
        and b["real_minus_shuffled_spearman"] is not None
        and b["target_spearman"] >= 0.30
        and b["real_minus_shuffled_spearman"] >= 0.20
        for b in blocks
    )

    rejected = all(
        b["target_spearman"] is not None
        and b["real_minus_shuffled_spearman"] is not None
        and b["target_spearman"] < 0.10
        and b["real_minus_shuffled_spearman"] < 0.10
        for b in blocks
    )

    if supported:
        decision = "UNIVERSAL_ALGEBRAIC_PRIOR_SUPPORTED"
    elif rejected:
        decision = "UNIVERSAL_ALGEBRAIC_PRIOR_REJECTED_FOR_THIS_FEATURE_CLASS"
    else:
        decision = "UNIVERSAL_ALGEBRAIC_PRIOR_UNRESOLVED"

    return {
        "version": "NABU_EXTERNAL_HIGHER_ORDER_PRIOR_TRANSFER_V1",
        "decision": decision,
        "root_cause_status": "ROOT_CAUSE_NOT_YET_IDENTIFIED",
        "gb1_to_trpb": {
            "spearman": gb_to_tr["target_spearman"],
            "real_minus_shuffled": gb_to_tr[
                "real_minus_shuffled_spearman"
            ],
        },
        "trpb_to_gb1": {
            "spearman": tr_to_gb["target_spearman"],
            "real_minus_shuffled": tr_to_gb[
                "real_minus_shuffled_spearman"
            ],
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

    trpb_path = Path(args.trpb)
    gb1_path = Path(args.gb1)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    trpb_data = s3.load_trpb(trpb_path)
    gb1_data = s3.load_gb1(gb1_path)

    trpb = build_triples(trpb_data)
    gb1 = build_triples(gb1_data)

    gb_to_tr = transfer(gb1, trpb)
    tr_to_gb = transfer(trpb, gb1)
    matrix = classify(gb_to_tr, tr_to_gb)

    (out / "GB1_TO_TRPB.json").write_text(
        json.dumps(gb_to_tr, indent=2), encoding="utf-8"
    )
    (out / "TRPB_TO_GB1.json").write_text(
        json.dumps(tr_to_gb, indent=2), encoding="utf-8"
    )
    (out / "PRIOR_TRANSFER_MATRIX.json").write_text(
        json.dumps(matrix, indent=2), encoding="utf-8"
    )

    manifest = {
        "version": "NABU_EXTERNAL_HIGHER_ORDER_PRIOR_TRANSFER_RUN_V1",
        "seed": SEED,
        "trpb_sha256": sha256_file(trpb_path),
        "gb1_sha256": sha256_file(gb1_path),
        "gb1_git_blob_sha": s3.s2.git_blob_sha(gb1_path),
        "nucb_consumed": False,
        "phase3_opened": False,
    }
    (out / "RUN_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )

    hashes = {}
    for name in (
        "GB1_TO_TRPB.json",
        "TRPB_TO_GB1.json",
        "PRIOR_TRANSFER_MATRIX.json",
        "RUN_MANIFEST.json",
    ):
        hashes[name] = sha256_file(out / name)
    (out / "OUTPUT_HASHES.json").write_text(
        json.dumps(hashes, indent=2), encoding="utf-8"
    )

    print(json.dumps(matrix, indent=2))


if __name__ == "__main__":
    main()
