from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from run_v9_pair_transfer_ired import derive_reference, load_ired, mutation_set
from nabu_protein.higher_order import fnv1a32


SEED = 161
N_FOLDS = 5
VERSION = "NABU_CLEAN_PAIR_UNSEEN_EDGE_POC_V1"


def assign_fold(candidate_id: str) -> int:
    return fnv1a32(candidate_id + "|NABU_CLEAN_PAIR_CV|161") % N_FOLDS


def safe_spearman(target: np.ndarray, prediction: np.ndarray) -> float | None:
    if len(target) < 2 or np.all(prediction == prediction[0]):
        return None
    value = float(spearmanr(target, prediction).statistic)
    return value if np.isfinite(value) else None


def safe_pearson(target: np.ndarray, prediction: np.ndarray) -> float | None:
    if len(target) < 2 or np.all(prediction == prediction[0]):
        return None
    value = float(pearsonr(target, prediction).statistic)
    return value if np.isfinite(value) else None


def evaluate(target: np.ndarray, prediction: np.ndarray) -> dict:
    target = np.asarray(target, dtype=float)
    prediction = np.asarray(prediction, dtype=float)
    if len(target) != len(prediction) or len(target) == 0:
        raise ValueError("Target/prediction length mismatch or empty evaluation set.")
    if not np.isfinite(target).all() or not np.isfinite(prediction).all():
        raise ValueError("Evaluation requires finite values.")

    error = prediction - target
    return {
        "count": int(len(target)),
        "spearman": safe_spearman(target, prediction),
        "pearson": safe_pearson(target, prediction),
        "rmse": float(np.sqrt(np.mean(error ** 2))),
        "mae": float(np.mean(np.abs(error))),
        "sign_accuracy": float(np.mean(np.sign(prediction) == np.sign(target))),
        "prediction_mean": float(np.mean(prediction)),
        "prediction_std": float(np.std(prediction)),
        "target_mean": float(np.mean(target)),
        "target_std": float(np.std(target)),
    }


def fit_node_additive(
    pairs: list[tuple[str, str]],
    residual: np.ndarray,
    *,
    permute: bool,
    fold: int,
) -> dict:
    residual = np.asarray(residual, dtype=float)
    if len(pairs) != len(residual) or len(pairs) == 0:
        raise ValueError("Pair/residual mismatch or empty training set.")

    nodes = sorted({mutation for pair in pairs for mutation in pair})
    node_index = {mutation: index for index, mutation in enumerate(nodes)}

    design = np.zeros((len(pairs), len(nodes) + 1), dtype=float)
    design[:, 0] = 1.0
    support = Counter()

    for row, pair in enumerate(pairs):
        for mutation in pair:
            design[row, 1 + node_index[mutation]] = 1.0
            support[mutation] += 1

    target = residual.copy()
    if permute:
        rng = np.random.default_rng(SEED + fold)
        target = rng.permutation(target)

    coefficient, _, rank, singular_values = np.linalg.lstsq(
        design,
        target,
        rcond=None,
    )

    return {
        "bias": float(coefficient[0]),
        "effect": {
            mutation: float(coefficient[1 + node_index[mutation]])
            for mutation in nodes
        },
        "support": {
            mutation: int(support[mutation])
            for mutation in nodes
        },
        "rank": int(rank),
        "parameter_count": int(len(nodes) + 1),
        "row_count": int(len(pairs)),
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
        "permuted": bool(permute),
    }


def predict_node_additive(
    pair: tuple[str, str],
    model: dict,
) -> float | None:
    left, right = pair
    effect = model["effect"]
    if left not in effect or right not in effect:
        return None
    return float(model["bias"] + effect[left] + effect[right])


def prediction_hash(frame: pd.DataFrame, column: str) -> str:
    payload = "\n".join(
        f"{row.candidate_id},{getattr(row, column):.17g}"
        for row in frame[["candidate_id", column]].itertuples(index=False)
        if np.isfinite(getattr(row, column))
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def build_clean_pair_table(source_gz: Path) -> tuple[dict, pd.DataFrame]:
    source_sha, fit_frame, validation_frame, test_frame = load_ired(source_gz)
    reference = derive_reference(fit_frame["sequence"].tolist())

    work = fit_frame[["sequence", "target"]].copy()
    work["mutation_set"] = [
        mutation_set(sequence, reference)
        for sequence in work["sequence"]
    ]
    work["mutation_count"] = work["mutation_set"].map(len)

    count_distribution = (
        work["mutation_count"]
        .value_counts()
        .sort_index()
        .to_dict()
    )

    wt_rows = work[work["mutation_count"] == 0]
    if len(wt_rows) != 1:
        raise RuntimeError(f"Expected exactly one WT row, found {len(wt_rows)}.")
    wt_target = float(wt_rows.iloc[0]["target"])

    singles = work[work["mutation_count"] == 1].copy()
    if singles["mutation_set"].duplicated().any():
        raise RuntimeError("Duplicate single-mutant identities found.")

    single_target = {
        mutations[0]: float(target)
        for mutations, target in zip(
            singles["mutation_set"],
            singles["target"],
        )
    }

    doubles = work[work["mutation_count"] == 2].copy()
    clean_rows = []
    missing_single_rows = []

    for row in doubles.itertuples(index=False):
        left, right = row.mutation_set
        missing = [
            mutation
            for mutation in (left, right)
            if mutation not in single_target
        ]
        record = {
            "candidate_id": str(row.sequence),
            "left": left,
            "right": right,
            "pair": (left, right),
            "double_target": float(row.target),
            "missing_single_count": int(len(missing)),
        }

        if missing:
            missing_single_rows.append(record)
            continue

        additive_expectation = (
            single_target[left]
            + single_target[right]
            - wt_target
        )
        residual = float(row.target - additive_expectation)
        clean_rows.append(
            {
                **record,
                "left_single_target": float(single_target[left]),
                "right_single_target": float(single_target[right]),
                "wt_target": wt_target,
                "clean_additive_expectation": float(additive_expectation),
                "clean_pair_residual": residual,
                "fold": int(assign_fold(str(row.sequence))),
            }
        )

    clean = pd.DataFrame(clean_rows)
    if clean.empty:
        raise RuntimeError("No clean double-mutant rows with both singleton counterparts.")

    pair_keys = clean[["left", "right"]].apply(tuple, axis=1)
    if pair_keys.duplicated().any():
        raise RuntimeError("Duplicate clean pair identities found.")

    metadata = {
        "version": VERSION,
        "source_sha256": source_sha,
        "reference_sha256": hashlib.sha256(reference.encode("utf-8")).hexdigest(),
        "reference_length": int(len(reference)),
        "fit_count": int(len(fit_frame)),
        "validation_rows_unused": int(len(validation_frame)),
        "test_rows_unused": int(len(test_frame)),
        "fit_mutation_count_distribution": {
            str(key): int(value)
            for key, value in sorted(count_distribution.items())
        },
        "wt_target": wt_target,
        "single_count": int(len(singles)),
        "double_count": int(len(doubles)),
        "clean_pair_count": int(len(clean)),
        "double_missing_single_counterpart_count": int(len(missing_single_rows)),
        "clean_pair_fraction_of_doubles": float(len(clean) / len(doubles)),
        "clean_pair_unique_nodes": int(
            len(set(clean["left"]).union(set(clean["right"])))
        ),
        "fold_counts": {
            str(fold): int((clean["fold"] == fold).sum())
            for fold in range(N_FOLDS)
        },
        "clean_residual_summary": {
            "mean": float(clean["clean_pair_residual"].mean()),
            "std": float(clean["clean_pair_residual"].std(ddof=0)),
            "min": float(clean["clean_pair_residual"].min()),
            "max": float(clean["clean_pair_residual"].max()),
            "median": float(clean["clean_pair_residual"].median()),
        },
    }
    return metadata, clean


def run(source_gz: Path, output_dir: Path) -> dict:
    metadata, clean = build_clean_pair_table(source_gz)

    prediction_rows = []
    fold_results = []

    for fold in range(N_FOLDS):
        train = clean[clean["fold"] != fold].copy()
        hold = clean[clean["fold"] == fold].copy()

        train_pairs = list(zip(train["left"], train["right"]))
        train_target = train["clean_pair_residual"].to_numpy(dtype=float)
        train_mean = float(np.mean(train_target))

        model = fit_node_additive(
            train_pairs,
            train_target,
            permute=False,
            fold=fold,
        )
        permuted = fit_node_additive(
            train_pairs,
            train_target,
            permute=True,
            fold=fold,
        )

        train_nodes = set(model["effect"])

        fold_prediction_rows = []
        for row in hold.itertuples(index=False):
            left_seen = row.left in train_nodes
            right_seen = row.right in train_nodes
            if left_seen and right_seen:
                coverage_class = "NODE_COVERED"
            elif left_seen or right_seen:
                coverage_class = "ONE_NODE_COLD"
            else:
                coverage_class = "TWO_NODE_COLD"

            relation_prediction = predict_node_additive(
                (row.left, row.right),
                model,
            )
            permuted_prediction = predict_node_additive(
                (row.left, row.right),
                permuted,
            )

            fold_prediction_rows.append(
                {
                    "candidate_id": row.candidate_id,
                    "left": row.left,
                    "right": row.right,
                    "fold": int(fold),
                    "coverage_class": coverage_class,
                    "true_clean_pair_residual": float(row.clean_pair_residual),
                    "ZERO_INTERACTION": 0.0,
                    "TRAIN_MEAN": train_mean,
                    "NODE_ADDITIVE_RELATION": (
                        float(relation_prediction)
                        if relation_prediction is not None
                        else np.nan
                    ),
                    "PERMUTED_NODE_ADDITIVE_CONTROL": (
                        float(permuted_prediction)
                        if permuted_prediction is not None
                        else np.nan
                    ),
                }
            )

        fold_frame = pd.DataFrame(fold_prediction_rows)
        prediction_rows.extend(fold_prediction_rows)

        covered = fold_frame[
            fold_frame["coverage_class"] == "NODE_COVERED"
        ].copy()
        if covered.empty:
            raise RuntimeError(f"Fold {fold} has no NODE_COVERED held-out edges.")

        target = covered["true_clean_pair_residual"].to_numpy(dtype=float)
        fold_metrics = {
            name: evaluate(
                target,
                covered[name].to_numpy(dtype=float),
            )
            for name in (
                "ZERO_INTERACTION",
                "TRAIN_MEAN",
                "NODE_ADDITIVE_RELATION",
                "PERMUTED_NODE_ADDITIVE_CONTROL",
            )
        }

        fold_results.append(
            {
                "fold": int(fold),
                "train_count": int(len(train)),
                "hold_count": int(len(hold)),
                "node_covered_count": int(len(covered)),
                "one_node_cold_count": int(
                    (fold_frame["coverage_class"] == "ONE_NODE_COLD").sum()
                ),
                "two_node_cold_count": int(
                    (fold_frame["coverage_class"] == "TWO_NODE_COLD").sum()
                ),
                "train_node_count": int(len(train_nodes)),
                "design_rank": int(model["rank"]),
                "parameter_count": int(model["parameter_count"]),
                "rank_deficiency": int(
                    model["parameter_count"] - model["rank"]
                ),
                "metrics": fold_metrics,
            }
        )

    predictions = pd.DataFrame(prediction_rows).sort_values(
        ["fold", "candidate_id"],
        kind="stable",
    ).reset_index(drop=True)

    covered = predictions[
        predictions["coverage_class"] == "NODE_COVERED"
    ].copy()
    target = covered["true_clean_pair_residual"].to_numpy(dtype=float)

    pooled = {
        name: evaluate(
            target,
            covered[name].to_numpy(dtype=float),
        )
        for name in (
            "ZERO_INTERACTION",
            "TRAIN_MEAN",
            "NODE_ADDITIVE_RELATION",
            "PERMUTED_NODE_ADDITIVE_CONTROL",
        )
    }

    relation = pooled["NODE_ADDITIVE_RELATION"]
    zero = pooled["ZERO_INTERACTION"]
    mean = pooled["TRAIN_MEAN"]
    permuted = pooled["PERMUTED_NODE_ADDITIVE_CONTROL"]

    fold_direction_wins = 0
    for fold_result in fold_results:
        relation_rho = fold_result["metrics"]["NODE_ADDITIVE_RELATION"]["spearman"]
        permuted_rho = fold_result["metrics"][
            "PERMUTED_NODE_ADDITIVE_CONTROL"
        ]["spearman"]
        if (
            relation_rho is not None
            and permuted_rho is not None
            and relation_rho > permuted_rho
        ):
            fold_direction_wins += 1

    checks = {
        "positive_spearman": bool(
            relation["spearman"] is not None
            and relation["spearman"] > 0.0
        ),
        "rmse_beats_zero_interaction": bool(
            relation["rmse"] < zero["rmse"]
        ),
        "rmse_beats_train_mean": bool(
            relation["rmse"] < mean["rmse"]
        ),
        "spearman_beats_permuted_control": bool(
            relation["spearman"] is not None
            and permuted["spearman"] is not None
            and relation["spearman"] > permuted["spearman"]
        ),
        "sign_accuracy_above_half": bool(
            relation["sign_accuracy"] > 0.5
        ),
        "fold_direction_consistency_at_least_4_of_5": bool(
            fold_direction_wins >= 4
        ),
    }

    result = {
        **metadata,
        "status": "POST_PHASE2_DEVELOPMENT_DIAGNOSTIC",
        "phase3_opened": False,
        "scientific_claim": False,
        "primary_population": "NODE_COVERED_UNSEEN_DOUBLE_EDGES",
        "pooled_node_covered_count": int(len(covered)),
        "pooled_one_node_cold_count": int(
            (predictions["coverage_class"] == "ONE_NODE_COLD").sum()
        ),
        "pooled_two_node_cold_count": int(
            (predictions["coverage_class"] == "TWO_NODE_COLD").sum()
        ),
        "pooled_metrics": pooled,
        "rmse_improvement_over_zero": float(
            zero["rmse"] - relation["rmse"]
        ),
        "rmse_improvement_over_train_mean": float(
            mean["rmse"] - relation["rmse"]
        ),
        "fold_direction_wins_vs_permuted": int(fold_direction_wins),
        "fold_results": fold_results,
        "prediction_sha256": {
            name: prediction_hash(covered, name)
            for name in (
                "ZERO_INTERACTION",
                "TRAIN_MEAN",
                "NODE_ADDITIVE_RELATION",
                "PERMUTED_NODE_ADDITIVE_CONTROL",
            )
        },
        "decision_checks": checks,
        "development_decision": (
            "IDENTITY_ONLY_CLEAN_PAIR_TRANSFER_SUPPORTED"
            if all(checks.values())
            else "IDENTITY_ONLY_CLEAN_PAIR_TRANSFER_NOT_SUPPORTED"
        ),
        "interpretation_boundary": (
            "This POC uses only the already-revealed IRED fit/train partition. "
            "It is diagnostic development evidence and does not alter Phase 2."
        ),
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    predictions.to_csv(
        output_dir / "CLEAN_PAIR_OOF_PREDICTIONS.csv",
        index=False,
    )
    (output_dir / "CLEAN_PAIR_POC_RESULTS.json").write_text(
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
