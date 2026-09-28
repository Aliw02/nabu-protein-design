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
ROOT = HERE.parent
FALS = ROOT / "root_cause_falsification"
PRIOR = ROOT / "root_cause_prior_transfer"
sys.path.insert(0, str(FALS))
sys.path.insert(0, str(PRIOR))

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
    value = float(spearmanr(a, b).statistic)
    return value if np.isfinite(value) else None


def finite_pearson(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) < 2 or np.std(a) <= EPS or np.std(b) <= EPS:
        return None
    value = float(pearsonr(a, b).statistic)
    return value if np.isfinite(value) else None


def build_order_table(dataset, order: int):
    transform = s1.RobustAsinhTransform()
    transform.fit(dataset["fit_y"])

    fit_y = transform.transform(dataset["fit_y"])
    test_y = transform.transform(dataset["test_y"])

    lookup, duplicates = s1.unique_lookup(
        dataset["fit_sets"], fit_y, f"{dataset['name']}_cross_order_fit"
    )
    if duplicates:
        raise RuntimeError(
            f"{dataset['name']}: duplicate low-order mutation sets."
        )

    wt, e1, e2 = s1.exact_components(lookup)
    sigma_low = float(np.std(fit_y))
    if sigma_low <= EPS:
        raise RuntimeError(f"{dataset['name']}: zero low-order target scale.")

    X = []
    residual = []
    o2_values = []
    transformed_target = []
    ids = []
    mutation_sets = []

    for cid, muts, target, observed_order in zip(
        dataset["test_ids"],
        dataset["test_sets"],
        test_y,
        dataset["orders"],
    ):
        if int(observed_order) != int(order):
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
        r = float(target - o2)

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
        residual.append(r / sigma_low)
        o2_values.append(o2)
        transformed_target.append(float(target))
        ids.append(str(cid))
        mutation_sets.append(tuple(muts))

    if not X:
        raise RuntimeError(
            f"{dataset['name']}: no exact-O2-supported order-{order} rows."
        )

    return {
        "dataset": dataset["name"],
        "order": int(order),
        "X": np.asarray(X, dtype=float),
        "residual": np.asarray(residual, dtype=float),
        "o2": np.asarray(o2_values, dtype=float),
        "target": np.asarray(transformed_target, dtype=float),
        "ids": np.asarray(ids, dtype=str),
        "mutation_sets": mutation_sets,
        "sigma_low": sigma_low,
        "transform": {
            "center": float(transform.center),
            "scale": float(transform.scale),
            "scale_source": transform.scale_source,
        },
    }


def fit_rf(X, y):
    model = RandomForestRegressor(
        n_estimators=300,
        min_samples_leaf=50,
        max_features="sqrt",
        random_state=SEED,
        n_jobs=1,
    )
    model.fit(X, y)
    return model


def fit_affine(o2_dimensionless, residual_dimensionless):
    x = np.asarray(o2_dimensionless, dtype=float)
    y = np.asarray(residual_dimensionless, dtype=float)
    A = np.column_stack([np.ones(len(x), dtype=float), x])
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    return float(beta[0]), float(beta[1])


def target_b2_reference(dataset, target_table):
    model = s1.NabuV83Model().fit(
        mutation_sets=dataset["fit_sets"],
        labels=dataset["fit_y"],
        candidate_ids=dataset["fit_ids"],
    )
    base = model.model["base"]
    global_mean = float(base["global_mean"])

    scores = []
    for muts in target_table["mutation_sets"]:
        if all(m in base["main"] for m in muts):
            scores.append(
                float(s1.additive_score(muts, global_mean, base["main"]))
            )
        else:
            scores.append(global_mean)

    return np.asarray(scores, dtype=float)


def transfer_cross_order(source_dataset, target_dataset):
    source_triples = build_order_table(source_dataset, 3)
    target_quads = build_order_table(target_dataset, 4)

    real_model = fit_rf(
        source_triples["X"],
        source_triples["residual"],
    )
    real_pred_residual = np.asarray(
        real_model.predict(target_quads["X"]),
        dtype=float,
    )

    rng = np.random.default_rng(SEED)
    shuffled_source_residual = rng.permutation(source_triples["residual"])
    shuffled_model = fit_rf(
        source_triples["X"],
        shuffled_source_residual,
    )
    shuffled_pred_residual = np.asarray(
        shuffled_model.predict(target_quads["X"]),
        dtype=float,
    )

    source_o2_dimless = (
        source_triples["o2"] / source_triples["sigma_low"]
    )
    target_o2_dimless = (
        target_quads["o2"] / target_quads["sigma_low"]
    )
    affine_intercept, affine_slope = fit_affine(
        source_o2_dimless,
        source_triples["residual"],
    )
    affine_pred_residual = (
        affine_intercept + affine_slope * target_o2_dimless
    )

    true_residual = target_quads["residual"]

    real_s = finite_spearman(true_residual, real_pred_residual)
    shuffled_s = finite_spearman(true_residual, shuffled_pred_residual)
    affine_s = finite_spearman(true_residual, affine_pred_residual)

    sigma_target = target_quads["sigma_low"]
    reconstructed_real = (
        target_quads["o2"] + real_pred_residual * sigma_target
    )
    reconstructed_shuffled = (
        target_quads["o2"] + shuffled_pred_residual * sigma_target
    )
    reconstructed_affine = (
        target_quads["o2"] + affine_pred_residual * sigma_target
    )

    target_true = target_quads["target"]
    b2 = target_b2_reference(target_dataset, target_quads)

    fitness_real_s = finite_spearman(target_true, reconstructed_real)
    fitness_o2_s = finite_spearman(target_true, target_quads["o2"])
    fitness_b2_s = finite_spearman(target_true, b2)
    fitness_shuffled_s = finite_spearman(
        target_true, reconstructed_shuffled
    )
    fitness_affine_s = finite_spearman(
        target_true, reconstructed_affine
    )

    return {
        "source": source_dataset["name"],
        "target": target_dataset["name"],
        "source_order": 3,
        "target_order": 4,
        "source_supported_triples": int(len(source_triples["residual"])),
        "target_supported_quads": int(len(target_quads["residual"])),
        "residual_prediction": {
            "spearman": real_s,
            "pearson": finite_pearson(
                true_residual, real_pred_residual
            ),
            "rmse": float(
                np.sqrt(
                    np.mean(
                        (true_residual - real_pred_residual) ** 2
                    )
                )
            ),
            "shuffled_control_spearman": shuffled_s,
            "source_affine_control_spearman": affine_s,
            "real_minus_shuffled_spearman": (
                None
                if real_s is None or shuffled_s is None
                else float(real_s - shuffled_s)
            ),
            "real_minus_source_affine_spearman": (
                None
                if real_s is None or affine_s is None
                else float(real_s - affine_s)
            ),
        },
        "reconstructed_fitness": {
            "real_prior_spearman": fitness_real_s,
            "exact_o2_spearman": fitness_o2_s,
            "target_b2_spearman": fitness_b2_s,
            "shuffled_prior_spearman": fitness_shuffled_s,
            "source_affine_spearman": fitness_affine_s,
            "real_minus_exact_o2": (
                None
                if fitness_real_s is None or fitness_o2_s is None
                else float(fitness_real_s - fitness_o2_s)
            ),
            "real_minus_target_b2": (
                None
                if fitness_real_s is None or fitness_b2_s is None
                else float(fitness_real_s - fitness_b2_s)
            ),
        },
        "source_affine_control": {
            "intercept": affine_intercept,
            "slope": affine_slope,
        },
        "source_sigma_low": float(source_triples["sigma_low"]),
        "target_sigma_low": float(target_quads["sigma_low"]),
        "source_transform": source_triples["transform"],
        "target_transform": target_quads["transform"],
    }


def replicate_triple_transfer(source_dataset, target_dataset):
    source_triples = build_order_table(source_dataset, 3)
    target_triples = build_order_table(target_dataset, 3)

    model = fit_rf(
        source_triples["X"],
        source_triples["residual"],
    )
    pred = np.asarray(model.predict(target_triples["X"]), dtype=float)

    rng = np.random.default_rng(SEED)
    shuffled = rng.permutation(source_triples["residual"])
    control = fit_rf(source_triples["X"], shuffled)
    control_pred = np.asarray(
        control.predict(target_triples["X"]),
        dtype=float,
    )

    real_s = finite_spearman(target_triples["residual"], pred)
    shuffled_s = finite_spearman(
        target_triples["residual"], control_pred
    )

    return {
        "source": source_dataset["name"],
        "target": target_dataset["name"],
        "source_supported_triples": int(len(source_triples["residual"])),
        "target_supported_triples": int(len(target_triples["residual"])),
        "residual_spearman": real_s,
        "shuffled_control_spearman": shuffled_s,
        "real_minus_shuffled_spearman": (
            None
            if real_s is None or shuffled_s is None
            else float(real_s - shuffled_s)
        ),
    }


def direction_pass(block):
    r = block["residual_prediction"]
    f = block["reconstructed_fitness"]

    values = (
        r["spearman"],
        r["real_minus_shuffled_spearman"],
        r["real_minus_source_affine_spearman"],
        f["real_prior_spearman"],
        f["exact_o2_spearman"],
        f["target_b2_spearman"],
    )
    if any(v is None for v in values):
        return False

    return bool(
        r["spearman"] >= 0.30
        and r["real_minus_shuffled_spearman"] >= 0.20
        and r["real_minus_source_affine_spearman"] >= 0.10
        and f["real_prior_spearman"] >= f["exact_o2_spearman"] + 0.10
        and f["real_prior_spearman"] >= f["target_b2_spearman"] - 0.02
    )


def direction_rejected(block):
    r = block["residual_prediction"]
    if (
        r["spearman"] is None
        or r["real_minus_shuffled_spearman"] is None
    ):
        return False
    return bool(
        r["spearman"] < 0.10
        and r["real_minus_shuffled_spearman"] < 0.10
    )


def classify(gb_to_tr, tr_to_gb, rep_gb_to_tr, rep_tr_to_gb):
    passes = {
        "GB1_triples_to_TrpB_quads": direction_pass(gb_to_tr),
        "TrpB_triples_to_GB1_quads": direction_pass(tr_to_gb),
    }

    if all(passes.values()):
        decision = "CROSS_LANDSCAPE_CROSS_ORDER_PRIOR_SUPPORTED"
    elif (
        direction_rejected(gb_to_tr)
        and direction_rejected(tr_to_gb)
    ):
        decision = (
            "CROSS_LANDSCAPE_CROSS_ORDER_PRIOR_"
            "REJECTED_FOR_THIS_FEATURE_CLASS"
        )
    else:
        decision = "CROSS_LANDSCAPE_CROSS_ORDER_PRIOR_UNRESOLVED"

    return {
        "version": "NABU_CROSS_LANDSCAPE_CROSS_ORDER_PRIOR_V1",
        "decision": decision,
        "root_cause_status": "ROOT_CAUSE_NOT_YET_IDENTIFIED",
        "direction_gate": passes,
        "GB1_triples_to_TrpB_quads": {
            "residual_spearman": gb_to_tr["residual_prediction"]["spearman"],
            "real_minus_shuffled": gb_to_tr[
                "residual_prediction"
            ]["real_minus_shuffled_spearman"],
            "real_minus_source_affine": gb_to_tr[
                "residual_prediction"
            ]["real_minus_source_affine_spearman"],
            "reconstructed_fitness_spearman": gb_to_tr[
                "reconstructed_fitness"
            ]["real_prior_spearman"],
            "exact_o2_spearman": gb_to_tr[
                "reconstructed_fitness"
            ]["exact_o2_spearman"],
            "target_b2_spearman": gb_to_tr[
                "reconstructed_fitness"
            ]["target_b2_spearman"],
        },
        "TrpB_triples_to_GB1_quads": {
            "residual_spearman": tr_to_gb["residual_prediction"]["spearman"],
            "real_minus_shuffled": tr_to_gb[
                "residual_prediction"
            ]["real_minus_shuffled_spearman"],
            "real_minus_source_affine": tr_to_gb[
                "residual_prediction"
            ]["real_minus_source_affine_spearman"],
            "reconstructed_fitness_spearman": tr_to_gb[
                "reconstructed_fitness"
            ]["real_prior_spearman"],
            "exact_o2_spearman": tr_to_gb[
                "reconstructed_fitness"
            ]["exact_o2_spearman"],
            "target_b2_spearman": tr_to_gb[
                "reconstructed_fitness"
            ]["target_b2_spearman"],
        },
        "triple_transfer_replication": {
            "GB1_to_TrpB": rep_gb_to_tr,
            "TrpB_to_GB1": rep_tr_to_gb,
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

    trpb = s3.load_trpb(trpb_path)
    gb1 = s3.load_gb1(gb1_path)

    gb_to_tr = transfer_cross_order(gb1, trpb)
    tr_to_gb = transfer_cross_order(trpb, gb1)

    rep_gb_to_tr = replicate_triple_transfer(gb1, trpb)
    rep_tr_to_gb = replicate_triple_transfer(trpb, gb1)

    matrix = classify(
        gb_to_tr,
        tr_to_gb,
        rep_gb_to_tr,
        rep_tr_to_gb,
    )

    (out / "GB1_TRIPLES_TO_TRPB_QUADS.json").write_text(
        json.dumps(gb_to_tr, indent=2),
        encoding="utf-8",
    )
    (out / "TRPB_TRIPLES_TO_GB1_QUADS.json").write_text(
        json.dumps(tr_to_gb, indent=2),
        encoding="utf-8",
    )
    (out / "CROSS_ORDER_MATRIX.json").write_text(
        json.dumps(matrix, indent=2),
        encoding="utf-8",
    )

    manifest = {
        "version": "NABU_CROSS_LANDSCAPE_CROSS_ORDER_PRIOR_RUN_V1",
        "seed": SEED,
        "trpb_sha256": sha256_file(trpb_path),
        "gb1_sha256": sha256_file(gb1_path),
        "gb1_git_blob_sha": s3.s2.git_blob_sha(gb1_path),
        "nucb_consumed": False,
        "phase3_opened": False,
    }
    (out / "RUN_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )

    hashes = {}
    for name in (
        "GB1_TRIPLES_TO_TRPB_QUADS.json",
        "TRPB_TRIPLES_TO_GB1_QUADS.json",
        "CROSS_ORDER_MATRIX.json",
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
