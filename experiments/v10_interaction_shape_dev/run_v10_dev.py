from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from itertools import combinations
from pathlib import Path

import joblib
import numpy as np
from scipy.stats import rankdata
from sklearn.ensemble import RandomForestRegressor

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
FALS = ROOT / "root_cause_falsification"
PAIR = ROOT / "v9_pair_transfer_dev"
sys.path.insert(0, str(FALS))
sys.path.insert(0, str(PAIR))

import run_falsification as s1
import run_stage3 as s3
from run_v9_pair_transfer_ired import evaluate
from nabu_protein.higher_order import additive_score
from nabu_protein.v83 import NabuV83Model

SEED = 161
EPS = 1e-12

FEATURE_NAMES = [
    "mean_abs_e1",
    "max_abs_e1",
    "std_e1",
    "min_e1",
    "max_e1",
    "mean_abs_e2",
    "max_abs_e2",
    "std_e2",
    "min_e2",
    "max_e2",
    "fraction_positive_e2",
]

RF_PARAMS = {
    "n_estimators": 300,
    "min_samples_leaf": 50,
    "max_features": "sqrt",
    "random_state": SEED,
    "n_jobs": 1,
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def percentile_rank(values):
    values = np.asarray(values, dtype=float)
    if len(values) == 1:
        return np.asarray([0.5], dtype=float)
    ranks = rankdata(values, method="average")
    return (ranks - 1.0) / float(len(values) - 1)


def fit_b2(dataset):
    model = NabuV83Model().fit(
        mutation_sets=dataset["fit_sets"],
        labels=dataset["fit_y"],
        candidate_ids=dataset["fit_ids"],
    )
    base = model.model["base"]
    global_mean = float(base["global_mean"])

    raw = np.empty(len(dataset["test_sets"]), dtype=float)
    support = np.zeros(len(dataset["test_sets"]), dtype=bool)

    for i, muts in enumerate(dataset["test_sets"]):
        ok = all(m in base["main"] for m in muts)
        support[i] = ok
        if ok:
            raw[i] = float(
                additive_score(muts, global_mean, base["main"])
            )
        else:
            raw[i] = global_mean

    return model, raw, support


def interaction_shape_features(dataset):
    transform = s1.RobustAsinhTransform()
    transform.fit(dataset["fit_y"])

    fit_y = transform.transform(dataset["fit_y"])
    lookup, duplicates = s1.unique_lookup(
        dataset["fit_sets"], fit_y, f"{dataset['name']}_v10_fit"
    )
    if duplicates:
        raise RuntimeError(
            f"{dataset['name']}: duplicate low-order mutation sets."
        )

    wt, e1, e2 = s1.exact_components(lookup)
    sigma_low = float(np.std(fit_y))
    if sigma_low <= EPS:
        raise RuntimeError(f"{dataset['name']}: zero low-order scale.")

    n = len(dataset["test_sets"])
    supported = np.zeros(n, dtype=bool)
    shape = np.full((n, len(FEATURE_NAMES)), np.nan, dtype=float)
    o2_only = np.full((n, 1), np.nan, dtype=float)

    for i, muts in enumerate(dataset["test_sets"]):
        if not all(m in e1 for m in muts):
            continue

        pairs = [tuple(pair) for pair in combinations(muts, 2)]
        if not all(pair in e2 for pair in pairs):
            continue

        e1v = np.asarray([e1[m] for m in muts], dtype=float)
        e2v = np.asarray([e2[p] for p in pairs], dtype=float)

        o1 = float(wt + np.sum(e1v))
        o2 = float(o1 + np.sum(e2v))

        supported[i] = True
        shape[i] = np.asarray(
            [
                float(np.mean(np.abs(e1v))) / sigma_low,
                float(np.max(np.abs(e1v))) / sigma_low,
                float(np.std(e1v)) / sigma_low,
                float(np.min(e1v)) / sigma_low,
                float(np.max(e1v)) / sigma_low,
                float(np.mean(np.abs(e2v))) / sigma_low,
                float(np.max(np.abs(e2v))) / sigma_low,
                float(np.std(e2v)) / sigma_low,
                float(np.min(e2v)) / sigma_low,
                float(np.max(e2v)) / sigma_low,
                float(np.mean(e2v > 0.0)),
            ],
            dtype=float,
        )
        o2_only[i, 0] = float(o2 / sigma_low)

    if not np.isfinite(shape[supported]).all():
        raise RuntimeError(f"{dataset['name']}: non-finite shape features.")
    if not np.isfinite(o2_only[supported]).all():
        raise RuntimeError(f"{dataset['name']}: non-finite O2-only features.")

    return {
        "shape": shape,
        "o2_only": o2_only,
        "supported": supported,
        "sigma_low": sigma_low,
        "transform": {
            "center": float(transform.center),
            "scale": float(transform.scale),
            "scale_source": transform.scale_source,
        },
        "component_counts": {
            "e1": int(len(e1)),
            "e2": int(len(e2)),
        },
    }


def prepare_landscape(dataset):
    _, b2_raw, b2_support = fit_b2(dataset)
    features = interaction_shape_features(dataset)

    b2_pct = percentile_rank(b2_raw)
    true_pct = percentile_rank(dataset["test_y"])
    rank_residual = true_pct - b2_pct

    return {
        "name": dataset["name"],
        "dataset": dataset,
        "b2_raw": b2_raw,
        "b2_percentile": b2_pct,
        "b2_support": b2_support,
        "true_percentile": true_pct,
        "rank_residual": rank_residual,
        "shape": features["shape"],
        "o2_only": features["o2_only"],
        "v10_supported": features["supported"],
        "sigma_low": features["sigma_low"],
        "transform": features["transform"],
        "component_counts": features["component_counts"],
    }


def fit_rf(X, y):
    model = RandomForestRegressor(**RF_PARAMS)
    model.fit(X, y)
    return model


def evaluate_region(target_y, prediction, mask):
    mask = np.asarray(mask, dtype=bool)
    if int(np.sum(mask)) < 2:
        return {"count": int(np.sum(mask))}
    return evaluate(
        np.asarray(target_y, dtype=float)[mask],
        np.asarray(prediction, dtype=float)[mask],
    )


def score_direction(source, target):
    source_mask = source["v10_supported"]
    target_mask = target["v10_supported"]

    if int(np.sum(source_mask)) < 100:
        raise RuntimeError(
            f"{source['name']}: insufficient V10-supported source rows."
        )
    if int(np.sum(target_mask)) < 100:
        raise RuntimeError(
            f"{target['name']}: insufficient V10-supported target rows."
        )

    source_y = source["rank_residual"][source_mask]

    shape_model = fit_rf(source["shape"][source_mask], source_y)
    shape_pred_supported = np.asarray(
        shape_model.predict(target["shape"][target_mask]),
        dtype=float,
    )

    shape_replay = fit_rf(source["shape"][source_mask], source_y)
    replay_pred = np.asarray(
        shape_replay.predict(target["shape"][target_mask]),
        dtype=float,
    )
    replay_max_abs = float(
        np.max(np.abs(shape_pred_supported - replay_pred))
    )

    o2_model = fit_rf(source["o2_only"][source_mask], source_y)
    o2_pred_supported = np.asarray(
        o2_model.predict(target["o2_only"][target_mask]),
        dtype=float,
    )

    shuffled_source_y = np.random.default_rng(SEED).permutation(source_y)
    shuffled_model = fit_rf(
        source["shape"][source_mask],
        shuffled_source_y,
    )
    shuffled_pred_supported = np.asarray(
        shuffled_model.predict(target["shape"][target_mask]),
        dtype=float,
    )

    shape_correction = np.zeros(len(target["b2_percentile"]), dtype=float)
    o2_correction = np.zeros(len(target["b2_percentile"]), dtype=float)
    shuffled_correction = np.zeros(
        len(target["b2_percentile"]),
        dtype=float,
    )

    shape_correction[target_mask] = shape_pred_supported
    o2_correction[target_mask] = o2_pred_supported
    shuffled_correction[target_mask] = shuffled_pred_supported

    v10 = target["b2_percentile"] + shape_correction
    o2_control = target["b2_percentile"] + o2_correction
    shuffled_control = target["b2_percentile"] + shuffled_correction
    b2_score = target["b2_percentile"]

    target_y = target["dataset"]["test_y"]
    orders = target["dataset"]["orders"]

    full = {
        "B2": evaluate(target_y, b2_score),
        "V10": evaluate(target_y, v10),
        "O2_ONLY_CONTROL": evaluate(target_y, o2_control),
        "SHUFFLED_CONTROL": evaluate(target_y, shuffled_control),
    }

    by_order = {}
    for order in sorted(set(orders.tolist())):
        mask = orders == order
        if int(np.sum(mask)) < 2:
            continue
        by_order[str(order)] = {
            "B2": evaluate_region(target_y, b2_score, mask),
            "V10": evaluate_region(target_y, v10, mask),
            "O2_ONLY_CONTROL": evaluate_region(
                target_y, o2_control, mask
            ),
            "SHUFFLED_CONTROL": evaluate_region(
                target_y, shuffled_control, mask
            ),
        }

    supported = {
        "count": int(np.sum(target_mask)),
        "fraction": float(np.mean(target_mask)),
        "B2": evaluate_region(target_y, b2_score, target_mask),
        "V10": evaluate_region(target_y, v10, target_mask),
    }
    unsupported_mask = ~target_mask
    unsupported = {
        "count": int(np.sum(unsupported_mask)),
        "fraction": float(np.mean(unsupported_mask)),
        "B2": evaluate_region(target_y, b2_score, unsupported_mask),
        "V10": evaluate_region(target_y, v10, unsupported_mask),
        "max_abs_score_difference": (
            0.0
            if int(np.sum(unsupported_mask)) == 0
            else float(
                np.max(
                    np.abs(
                        v10[unsupported_mask]
                        - b2_score[unsupported_mask]
                    )
                )
            )
        ),
    }

    required_orders = ("3", "4")
    missing_orders = [
        order for order in required_orders if order not in by_order
    ]
    if missing_orders:
        raise RuntimeError(
            f"{target['name']}: missing required orders {missing_orders}"
        )

    gate = {
        "whole_spearman_strictly_beats_b2": bool(
            full["V10"]["spearman"] > full["B2"]["spearman"]
        ),
        "whole_ndcg_not_below_b2": bool(
            full["V10"]["ndcg"] >= full["B2"]["ndcg"]
        ),
        "whole_top1_hits_not_below_b2": bool(
            full["V10"]["top1_percent_hits"]
            >= full["B2"]["top1_percent_hits"]
        ),
        "order3_spearman_not_below_b2": bool(
            by_order["3"]["V10"]["spearman"]
            >= by_order["3"]["B2"]["spearman"]
        ),
        "order4_spearman_not_below_b2": bool(
            by_order["4"]["V10"]["spearman"]
            >= by_order["4"]["B2"]["spearman"]
        ),
        "whole_spearman_beats_o2_only_control": bool(
            full["V10"]["spearman"]
            > full["O2_ONLY_CONTROL"]["spearman"]
        ),
        "whole_spearman_beats_shuffled_control": bool(
            full["V10"]["spearman"]
            > full["SHUFFLED_CONTROL"]["spearman"]
        ),
        "unsupported_exactly_b2": bool(
            unsupported["max_abs_score_difference"] == 0.0
        ),
        "deterministic_replay_exact": bool(replay_max_abs == 0.0),
    }
    gate["all_required"] = bool(all(gate.values()))

    return {
        "source": source["name"],
        "target": target["name"],
        "source_training_count": int(np.sum(source_mask)),
        "target_supported_count": int(np.sum(target_mask)),
        "target_total_count": int(len(target_y)),
        "target_support_fraction": float(np.mean(target_mask)),
        "full": full,
        "by_order": by_order,
        "supported_region": supported,
        "unsupported_region": unsupported,
        "determinism": {
            "replay_max_abs_difference": replay_max_abs,
        },
        "gate": gate,
        "target_transform": target["transform"],
        "source_transform": source["transform"],
    }


def train_final_prior(prepared_landscapes, out_dir: Path):
    X_parts = []
    y_parts = []
    per_dataset = {}

    for prepared in prepared_landscapes:
        mask = prepared["v10_supported"]
        X_parts.append(prepared["shape"][mask])
        y_parts.append(prepared["rank_residual"][mask])
        per_dataset[prepared["name"]] = int(np.sum(mask))

    X = np.vstack(X_parts)
    y = np.concatenate(y_parts)

    model = fit_rf(X, y)
    model_path = out_dir / "FINAL_PRIOR_MODEL.joblib"
    joblib.dump(model, model_path)

    freeze = {
        "version": "NABU_V10_FINAL_EXTERNAL_PRIOR_FREEZE_V1",
        "feature_names": FEATURE_NAMES,
        "rf_params": RF_PARAMS,
        "training_rows_total": int(len(y)),
        "training_rows_by_dataset": per_dataset,
        "training_target": (
            "true_global_percentile_rank_minus_B2_global_percentile_rank"
        ),
        "unsupported_candidate_correction": 0.0,
        "score_rule": (
            "B2_global_percentile_rank_plus_predicted_rank_correction"
        ),
        "model_sha256": sha256_file(model_path),
        "nucb_consumed": False,
        "phase3_opened": False,
    }
    (out_dir / "FINAL_PRIOR_FREEZE.json").write_text(
        json.dumps(freeze, indent=2),
        encoding="utf-8",
    )
    return freeze


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

    trpb = prepare_landscape(trpb_data)
    gb1 = prepare_landscape(gb1_data)

    gb_to_trpb = score_direction(gb1, trpb)
    trpb_to_gb1 = score_direction(trpb, gb1)

    pass_both = bool(
        gb_to_trpb["gate"]["all_required"]
        and trpb_to_gb1["gate"]["all_required"]
    )

    decision = "V10_DEV_PASS" if pass_both else "V10_DEV_REJECT"

    matrix = {
        "version": "NABU_V10_INTERACTION_SHAPE_DEV_V1",
        "decision": decision,
        "architecture": "B2_PLUS_EXTERNAL_INTERACTION_SHAPE_RANK_CORRECTION",
        "direction_gate": {
            "GB1_to_TrpB": bool(gb_to_trpb["gate"]["all_required"]),
            "TrpB_to_GB1": bool(trpb_to_gb1["gate"]["all_required"]),
        },
        "kill_rule_applies": bool(not pass_both),
        "nucb_consumed": False,
        "phase3_opened": False,
    }

    (out / "GB1_TO_TRPB_V10.json").write_text(
        json.dumps(gb_to_trpb, indent=2),
        encoding="utf-8",
    )
    (out / "TRPB_TO_GB1_V10.json").write_text(
        json.dumps(trpb_to_gb1, indent=2),
        encoding="utf-8",
    )

    if pass_both:
        freeze = train_final_prior([gb1, trpb], out)
        matrix["final_prior_freeze"] = {
            "training_rows_total": freeze["training_rows_total"],
            "model_sha256": freeze["model_sha256"],
        }

    (out / "V10_DEV_MATRIX.json").write_text(
        json.dumps(matrix, indent=2),
        encoding="utf-8",
    )

    manifest = {
        "version": "NABU_V10_INTERACTION_SHAPE_DEV_RUN_V1",
        "seed": SEED,
        "feature_names": FEATURE_NAMES,
        "rf_params": RF_PARAMS,
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
    for path in sorted(out.iterdir()):
        if path.is_file() and path.name != "OUTPUT_HASHES.json":
            hashes[path.name] = sha256_file(path)
    (out / "OUTPUT_HASHES.json").write_text(
        json.dumps(hashes, indent=2),
        encoding="utf-8",
    )

    print(json.dumps(matrix, indent=2))


if __name__ == "__main__":
    main()
