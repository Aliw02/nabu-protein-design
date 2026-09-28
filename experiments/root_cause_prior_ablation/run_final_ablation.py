from __future__ import annotations

import argparse
import hashlib
import json
import sys
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

import run_stage3 as s3
import run_prior_transfer as prior

SEED = 161
EPS = 1e-12

FEATURE_ARMS = {
    "O2_ONLY": [1],
    "SHAPE_ONLY": [3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14],
    "FULL": list(range(15)),
}


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


def fit_model(X, y):
    model = RandomForestRegressor(
        n_estimators=300,
        min_samples_leaf=50,
        max_features="sqrt",
        random_state=SEED,
        n_jobs=1,
    )
    model.fit(X, y)
    return model


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


def evaluate_arm(source, target, indices):
    X_source = source["X"][:, indices]
    X_target = target["X"][:, indices]

    model = fit_model(X_source, source["y"])
    prediction = np.asarray(model.predict(X_target), dtype=float)

    true_conditional, true_affine = residualize_against_o2(
        target["o2"], target["y"]
    )
    pred_conditional, pred_affine = residualize_against_o2(
        target["o2"], prediction
    )

    return {
        "feature_count": int(len(indices)),
        "target_spearman": finite_spearman(target["y"], prediction),
        "target_pearson": finite_pearson(target["y"], prediction),
        "target_rmse": float(
            np.sqrt(np.mean((target["y"] - prediction) ** 2))
        ),
        "conditional_residual_spearman": finite_spearman(
            true_conditional, pred_conditional
        ),
        "true_e3_vs_o2_affine": true_affine,
        "prediction_vs_o2_affine": pred_affine,
    }


def transfer_direction(source, target):
    arms = {
        name: evaluate_arm(source, target, indices)
        for name, indices in FEATURE_ARMS.items()
    }

    full_cond = arms["FULL"]["conditional_residual_spearman"]
    o2_cond = arms["O2_ONLY"]["conditional_residual_spearman"]
    shape_target = arms["SHAPE_ONLY"]["target_spearman"]

    gate = {
        "full_conditional_at_least_0_30": bool(
            full_cond is not None and full_cond >= 0.30
        ),
        "full_minus_o2_conditional_at_least_0_10": bool(
            full_cond is not None
            and o2_cond is not None
            and full_cond - o2_cond >= 0.10
        ),
        "shape_only_target_at_least_0_20": bool(
            shape_target is not None and shape_target >= 0.20
        ),
    }
    gate["all_required"] = bool(all(gate.values()))

    return {
        "source": source["dataset"],
        "target": target["dataset"],
        "source_count": int(len(source["y"])),
        "target_count": int(len(target["y"])),
        "arms": arms,
        "gate": gate,
        "derived_differences": {
            "full_minus_o2_conditional": (
                None
                if full_cond is None or o2_cond is None
                else float(full_cond - o2_cond)
            ),
            "shape_only_target_spearman": shape_target,
            "full_conditional_residual_spearman": full_cond,
            "o2_only_conditional_residual_spearman": o2_cond,
        },
    }


def classify(gb_to_trpb, trpb_to_gb):
    both_pass = bool(
        gb_to_trpb["gate"]["all_required"]
        and trpb_to_gb["gate"]["all_required"]
    )

    if both_pass:
        decision = "TRANSFERABLE_INTERACTION_SHAPE_SUPPORTED"
        stopping_action = (
            "STOP_ROOT_CAUSE_DIAGNOSTICS_BUILD_ONE_DERIVED_ARCHITECTURE"
        )
        algebraic_prior_path = "KEEP"
        next_prior_class = None
    else:
        decision = "TRANSFERABLE_INTERACTION_SHAPE_NOT_SUPPORTED"
        stopping_action = (
            "STOP_ALGEBRAIC_PRIOR_PATH_MOVE_TO_SEQUENCE_AWARE_PRETRAINED_CONTEXT"
        )
        algebraic_prior_path = "KILL"
        next_prior_class = "SEQUENCE_AWARE_PRETRAINED_CONTEXT"

    return {
        "version": "NABU_PRIOR_TRANSFER_FINAL_ABLATION_V1",
        "decision": decision,
        "stopping_action": stopping_action,
        "algebraic_prior_path": algebraic_prior_path,
        "next_prior_class": next_prior_class,
        "direction_gate": {
            "GB1_to_TrpB": bool(gb_to_trpb["gate"]["all_required"]),
            "TrpB_to_GB1": bool(trpb_to_gb["gate"]["all_required"]),
        },
        "root_cause_status": "ROOT_CAUSE_NOT_YET_IDENTIFIED",
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

    trpb = prior.build_triples(trpb_data)
    gb1 = prior.build_triples(gb1_data)

    gb_to_trpb = transfer_direction(gb1, trpb)
    trpb_to_gb = transfer_direction(trpb, gb1)
    matrix = classify(gb_to_trpb, trpb_to_gb)

    (out / "GB1_TO_TRPB_ABLATION.json").write_text(
        json.dumps(gb_to_trpb, indent=2), encoding="utf-8"
    )
    (out / "TRPB_TO_GB1_ABLATION.json").write_text(
        json.dumps(trpb_to_gb, indent=2), encoding="utf-8"
    )
    (out / "PRIOR_ABLATION_MATRIX.json").write_text(
        json.dumps(matrix, indent=2), encoding="utf-8"
    )

    manifest = {
        "version": "NABU_PRIOR_TRANSFER_FINAL_ABLATION_RUN_V1",
        "seed": SEED,
        "feature_arms": FEATURE_ARMS,
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
        "GB1_TO_TRPB_ABLATION.json",
        "TRPB_TO_GB1_ABLATION.json",
        "PRIOR_ABLATION_MATRIX.json",
        "RUN_MANIFEST.json",
    ):
        hashes[name] = sha256_file(out / name)

    (out / "OUTPUT_HASHES.json").write_text(
        json.dumps(hashes, indent=2), encoding="utf-8"
    )

    print(json.dumps(matrix, indent=2))


if __name__ == "__main__":
    main()
