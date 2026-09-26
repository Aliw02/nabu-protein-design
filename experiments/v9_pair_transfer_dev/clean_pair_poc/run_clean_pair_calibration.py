from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from run_clean_pair_poc import (
    N_FOLDS,
    SEED,
    build_clean_pair_table,
    evaluate,
    fit_node_additive,
    predict_node_additive,
)
from nabu_protein.higher_order import fnv1a32


VERSION = "NABU_CLEAN_PAIR_NESTED_OOF_CALIBRATION_V1"


def inner_fold(candidate_id: str, outer_fold: int) -> int:
    token = f"{candidate_id}|NABU_CLEAN_PAIR_INNER_CAL|{SEED}|{outer_fold}"
    return fnv1a32(token) % N_FOLDS


def fit_affine(raw_prediction: np.ndarray, target: np.ndarray) -> dict:
    raw_prediction = np.asarray(raw_prediction, dtype=float)
    target = np.asarray(target, dtype=float)
    if len(raw_prediction) != len(target) or len(target) < 2:
        raise ValueError("Calibration requires matched arrays with at least two rows.")
    if not np.isfinite(raw_prediction).all() or not np.isfinite(target).all():
        raise ValueError("Calibration requires finite arrays.")

    design = np.column_stack(
        [
            np.ones(len(raw_prediction), dtype=float),
            raw_prediction,
        ]
    )
    coefficient, _, rank, singular_values = np.linalg.lstsq(
        design,
        target,
        rcond=None,
    )
    return {
        "alpha": float(coefficient[0]),
        "beta": float(coefficient[1]),
        "rank": int(rank),
        "row_count": int(len(target)),
        "singular_min": (
            float(np.min(singular_values))
            if len(singular_values)
            else 0.0
        ),
        "singular_max": (
            float(np.max(singular_values))
            if len(singular_values)
            else 0.0
        ),
    }


def generate_inner_oof_calibration(
    outer_train: pd.DataFrame,
    outer_fold: int,
) -> tuple[dict, pd.DataFrame]:
    inner = outer_train.copy()
    inner["inner_fold"] = [
        inner_fold(candidate_id, outer_fold)
        for candidate_id in inner["candidate_id"]
    ]

    rows = []
    for fold in range(N_FOLDS):
        fit_frame = inner[inner["inner_fold"] != fold]
        hold_frame = inner[inner["inner_fold"] == fold]

        train_pairs = list(zip(fit_frame["left"], fit_frame["right"]))
        train_target = fit_frame["clean_pair_residual"].to_numpy(dtype=float)
        model = fit_node_additive(
            train_pairs,
            train_target,
            permute=False,
            fold=outer_fold * N_FOLDS + fold,
        )
        train_nodes = set(model["effect"])

        for row in hold_frame.itertuples(index=False):
            if row.left not in train_nodes or row.right not in train_nodes:
                continue
            prediction = predict_node_additive(
                (row.left, row.right),
                model,
            )
            if prediction is None:
                continue
            rows.append(
                {
                    "candidate_id": row.candidate_id,
                    "outer_fold": int(outer_fold),
                    "inner_fold": int(fold),
                    "raw_prediction": float(prediction),
                    "target": float(row.clean_pair_residual),
                }
            )

    oof = pd.DataFrame(rows)
    if len(oof) < 2:
        raise RuntimeError(
            f"Outer fold {outer_fold} produced insufficient inner OOF calibration rows."
        )

    calibration = fit_affine(
        oof["raw_prediction"].to_numpy(dtype=float),
        oof["target"].to_numpy(dtype=float),
    )
    return calibration, oof


def prediction_hash(frame: pd.DataFrame, column: str) -> str:
    payload = "\n".join(
        f"{row.candidate_id},{getattr(row, column):.17g}"
        for row in frame[["candidate_id", column]].itertuples(index=False)
        if np.isfinite(getattr(row, column))
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def run(source_gz: Path, output_dir: Path) -> dict:
    metadata, clean = build_clean_pair_table(source_gz)

    prediction_rows = []
    fold_results = []
    calibration_rows = []

    for outer_fold in range(N_FOLDS):
        outer_train = clean[clean["fold"] != outer_fold].copy()
        outer_hold = clean[clean["fold"] == outer_fold].copy()

        calibration, inner_oof = generate_inner_oof_calibration(
            outer_train,
            outer_fold,
        )

        train_pairs = list(zip(outer_train["left"], outer_train["right"]))
        train_target = outer_train["clean_pair_residual"].to_numpy(dtype=float)
        train_mean = float(np.mean(train_target))

        raw_model = fit_node_additive(
            train_pairs,
            train_target,
            permute=False,
            fold=outer_fold,
        )
        permuted_model = fit_node_additive(
            train_pairs,
            train_target,
            permute=True,
            fold=outer_fold,
        )
        train_nodes = set(raw_model["effect"])

        fold_rows = []
        for row in outer_hold.itertuples(index=False):
            node_covered = (
                row.left in train_nodes
                and row.right in train_nodes
            )
            if not node_covered:
                continue

            raw = predict_node_additive(
                (row.left, row.right),
                raw_model,
            )
            permuted = predict_node_additive(
                (row.left, row.right),
                permuted_model,
            )
            if raw is None or permuted is None:
                raise RuntimeError("NODE_COVERED row lacked a prediction.")

            calibrated = (
                calibration["alpha"]
                + calibration["beta"] * float(raw)
            )

            fold_rows.append(
                {
                    "candidate_id": row.candidate_id,
                    "outer_fold": int(outer_fold),
                    "true_clean_pair_residual": float(row.clean_pair_residual),
                    "ZERO_INTERACTION": 0.0,
                    "TRAIN_MEAN": train_mean,
                    "RAW_NODE_ADDITIVE_RELATION": float(raw),
                    "OOF_AFFINE_CALIBRATED_RELATION": float(calibrated),
                    "PERMUTED_RELATION_CONTROL": float(permuted),
                }
            )

        fold_frame = pd.DataFrame(fold_rows)
        if fold_frame.empty:
            raise RuntimeError(
                f"Outer fold {outer_fold} had no NODE_COVERED evaluation rows."
            )

        target = fold_frame["true_clean_pair_residual"].to_numpy(dtype=float)
        metrics = {
            name: evaluate(
                target,
                fold_frame[name].to_numpy(dtype=float),
            )
            for name in (
                "ZERO_INTERACTION",
                "TRAIN_MEAN",
                "RAW_NODE_ADDITIVE_RELATION",
                "OOF_AFFINE_CALIBRATED_RELATION",
                "PERMUTED_RELATION_CONTROL",
            )
        }

        prediction_rows.extend(fold_rows)
        calibration_rows.append(
            {
                "outer_fold": int(outer_fold),
                **calibration,
                "inner_oof_prediction_mean": float(
                    inner_oof["raw_prediction"].mean()
                ),
                "inner_oof_prediction_std": float(
                    inner_oof["raw_prediction"].std(ddof=0)
                ),
                "inner_oof_target_mean": float(inner_oof["target"].mean()),
                "inner_oof_target_std": float(
                    inner_oof["target"].std(ddof=0)
                ),
            }
        )
        fold_results.append(
            {
                "outer_fold": int(outer_fold),
                "outer_train_count": int(len(outer_train)),
                "outer_hold_count": int(len(outer_hold)),
                "node_covered_eval_count": int(len(fold_frame)),
                "calibration": calibration_rows[-1],
                "metrics": metrics,
            }
        )

    predictions = pd.DataFrame(prediction_rows).sort_values(
        ["outer_fold", "candidate_id"],
        kind="stable",
    ).reset_index(drop=True)

    target = predictions["true_clean_pair_residual"].to_numpy(dtype=float)
    arm_names = (
        "ZERO_INTERACTION",
        "TRAIN_MEAN",
        "RAW_NODE_ADDITIVE_RELATION",
        "OOF_AFFINE_CALIBRATED_RELATION",
        "PERMUTED_RELATION_CONTROL",
    )
    pooled = {
        name: evaluate(
            target,
            predictions[name].to_numpy(dtype=float),
        )
        for name in arm_names
    }

    calibrated = pooled["OOF_AFFINE_CALIBRATED_RELATION"]
    raw = pooled["RAW_NODE_ADDITIVE_RELATION"]
    zero = pooled["ZERO_INTERACTION"]
    train_mean = pooled["TRAIN_MEAN"]

    beta_all_positive = all(
        row["beta"] > 0.0
        for row in calibration_rows
    )
    fold_rmse_wins_vs_train_mean = sum(
        fold["metrics"]["OOF_AFFINE_CALIBRATED_RELATION"]["rmse"]
        < fold["metrics"]["TRAIN_MEAN"]["rmse"]
        for fold in fold_results
    )

    checks = {
        "calibrated_positive_spearman": bool(
            calibrated["spearman"] is not None
            and calibrated["spearman"] > 0.0
        ),
        "spearman_preserved_within_0_01": bool(
            calibrated["spearman"] is not None
            and raw["spearman"] is not None
            and calibrated["spearman"] >= raw["spearman"] - 0.01
        ),
        "calibrated_rmse_beats_zero": bool(
            calibrated["rmse"] < zero["rmse"]
        ),
        "calibrated_rmse_beats_train_mean": bool(
            calibrated["rmse"] < train_mean["rmse"]
        ),
        "calibrated_rmse_beats_raw": bool(
            calibrated["rmse"] < raw["rmse"]
        ),
        "calibrated_sign_accuracy_above_half": bool(
            calibrated["sign_accuracy"] > 0.5
        ),
        "beta_positive_all_folds": bool(beta_all_positive),
        "fold_rmse_beats_train_mean_at_least_4_of_5": bool(
            fold_rmse_wins_vs_train_mean >= 4
        ),
    }

    result = {
        "version": VERSION,
        "status": "POST_PHASE2_DEVELOPMENT_DIAGNOSTIC",
        "phase3_opened": False,
        "scientific_claim": False,
        "source_sha256": metadata["source_sha256"],
        "fit_count": metadata["fit_count"],
        "validation_rows_unused": metadata["validation_rows_unused"],
        "test_rows_unused": metadata["test_rows_unused"],
        "clean_pair_count": metadata["clean_pair_count"],
        "primary_population": "NODE_COVERED_UNSEEN_DOUBLE_EDGES",
        "pooled_count": int(len(predictions)),
        "single_calibration_mechanism": (
            "Nested-OOF ordinary least-squares affine calibration: "
            "epsilon = alpha + beta * raw_node_additive_prediction."
        ),
        "pooled_metrics": pooled,
        "calibration_by_outer_fold": calibration_rows,
        "fold_results": fold_results,
        "fold_rmse_wins_vs_train_mean": int(
            fold_rmse_wins_vs_train_mean
        ),
        "prediction_sha256": {
            name: prediction_hash(predictions, name)
            for name in arm_names
        },
        "decision_checks": checks,
        "development_decision": (
            "CLEAN_PAIR_AFFINE_CALIBRATION_SUPPORTED"
            if all(checks.values())
            else "CLEAN_PAIR_AFFINE_CALIBRATION_NOT_SUPPORTED"
        ),
        "next_step_if_fail": (
            "Close identity-only node-additive pair representation and "
            "move to richer relational sequence/structure features."
        ),
        "next_step_if_pass": (
            "Retain calibrated clean-pair relation as a candidate V9 component; "
            "do not open Phase 3 yet."
        ),
        "interpretation_boundary": (
            "Only the already-revealed IRED fit/train partition is used. "
            "Validation/test remain unused. This does not alter Phase 2."
        ),
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    predictions.to_csv(
        output_dir / "CLEAN_PAIR_CALIBRATION_OOF_PREDICTIONS.csv",
        index=False,
    )
    (output_dir / "CLEAN_PAIR_CALIBRATION_RESULTS.json").write_text(
        json.dumps(result, indent=2),
        encoding="utf-8",
    )

    print(json.dumps(result, indent=2))
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source_gz")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    run(Path(args.source_gz), Path(args.out))


if __name__ == "__main__":
    main()
