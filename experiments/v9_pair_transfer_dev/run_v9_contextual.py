from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

HERE = Path(__file__).resolve().parent
CLEAN_DIR = HERE / "clean_pair_poc"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(CLEAN_DIR))

from nabu_protein.v83 import NabuV83Model
from run_v9_pair_transfer_ired import (
    b2_predictions,
    derive_reference,
    evaluate,
    evaluate_by_mutation_count,
    fit_pair_factorization,
    load_ired,
    mutation_set,
    phase2_abstention_predictions,
)
from run_v9_1_degree_normalized import rank_hash, score_v91
from run_clean_pair_poc import build_clean_pair_table
from run_v9_candidate import (
    build_singleton_map,
    exact_pair_map,
    fit_global_oof_calibration,
    score_v9_clean_candidate,
)
from run_clean_pair_poc import fit_node_additive


VERSION = "NABU_V9_CONTEXTUAL_ESM_RESIDUAL_DEV_V1"
SEED = 161
MODEL_NAME = "facebook/esm2_t6_8M_UR50D"
MODEL_REVISION = "c731040fcd8d73dceaa04b0a8e6329b345b0f5df"
RIDGE_ALPHAS = np.asarray(
    [1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1000.0, 10000.0],
    dtype=float,
)
FEATURE_DIM = 640


def finite_spearman(target, prediction):
    value = float(spearmanr(target, prediction).statistic)
    return value if np.isfinite(value) else None


def sha256_array(values: np.ndarray, decimals: int = 8) -> str:
    array = np.asarray(values, dtype=np.float64)
    rounded = np.round(array, decimals=decimals)
    return hashlib.sha256(rounded.tobytes(order="C")).hexdigest()


def summarize_scale(values):
    values = np.asarray(values, dtype=float)
    return {
        "mean": float(np.mean(values)),
        "std": float(np.std(values)),
        "min": float(np.min(values)),
        "max": float(np.max(values)),
    }


def build_residual_training_table(
    fit_frame: pd.DataFrame,
    reference: str,
    clean_pairs: pd.DataFrame,
) -> pd.DataFrame:
    residual_by_sequence = {
        str(row.candidate_id): float(row.clean_pair_residual)
        for row in clean_pairs.itertuples(index=False)
    }

    rows = []
    excluded_double_count = 0

    for sequence, target in zip(
        fit_frame["sequence"].astype(str),
        fit_frame["target"].to_numpy(dtype=float),
    ):
        mutations = mutation_set(sequence, reference)
        count = len(mutations)

        if count == 0:
            rows.append(
                {
                    "sequence": sequence,
                    "mutation_count": 0,
                    "interaction_residual": 0.0,
                    "evidence_class": "WT",
                }
            )
        elif count == 1:
            rows.append(
                {
                    "sequence": sequence,
                    "mutation_count": 1,
                    "interaction_residual": 0.0,
                    "evidence_class": "SINGLE",
                }
            )
        elif count == 2 and sequence in residual_by_sequence:
            rows.append(
                {
                    "sequence": sequence,
                    "mutation_count": 2,
                    "interaction_residual": residual_by_sequence[sequence],
                    "evidence_class": "CLEAN_DOUBLE",
                }
            )
        elif count == 2:
            excluded_double_count += 1
        else:
            raise RuntimeError(
                "Unexpected mutation order >2 in IRED fit/train partition."
            )

    result = pd.DataFrame(rows)
    result.attrs["excluded_double_count"] = int(excluded_double_count)
    return result


def make_ridge_pipeline() -> Pipeline:
    inner_cv = KFold(
        n_splits=5,
        shuffle=True,
        random_state=SEED,
    )
    return Pipeline(
        [
            ("scale", StandardScaler()),
            (
                "ridge",
                RidgeCV(
                    alphas=RIDGE_ALPHAS,
                    cv=inner_cv,
                    scoring="neg_mean_squared_error",
                ),
            ),
        ]
    )


def fit_context_readout(features: np.ndarray, target: np.ndarray):
    features = np.asarray(features, dtype=float)
    target = np.asarray(target, dtype=float)

    if features.ndim != 2 or features.shape[1] != FEATURE_DIM:
        raise ValueError(
            f"Expected feature matrix (*,{FEATURE_DIM}), got {features.shape}."
        )
    if len(features) != len(target):
        raise ValueError("Feature/target row mismatch.")
    if not np.isfinite(features).all() or not np.isfinite(target).all():
        raise ValueError("Context readout requires finite inputs.")

    outer_cv = KFold(
        n_splits=5,
        shuffle=True,
        random_state=SEED + 101,
    )
    oof = cross_val_predict(
        make_ridge_pipeline(),
        features,
        target,
        cv=outer_cv,
        method="predict",
        n_jobs=1,
    )

    model = make_ridge_pipeline()
    model.fit(features, target)

    return model, np.asarray(oof, dtype=float)


def contextual_features_from_hidden(
    residue_hidden: np.ndarray,
    reference_hidden: np.ndarray,
    mutation_mask: np.ndarray,
) -> np.ndarray:
    residue_hidden = np.asarray(residue_hidden, dtype=np.float32)
    reference_hidden = np.asarray(reference_hidden, dtype=np.float32)
    mutation_mask = np.asarray(mutation_mask, dtype=np.float32)

    if residue_hidden.ndim != 3:
        raise ValueError("residue_hidden must be [batch, length, hidden].")
    if reference_hidden.shape != residue_hidden.shape[1:]:
        raise ValueError("Reference hidden shape mismatch.")
    if mutation_mask.shape != residue_hidden.shape[:2]:
        raise ValueError("Mutation mask shape mismatch.")

    delta = residue_hidden - reference_hidden[None, :, :]
    global_delta = delta.mean(axis=1)

    count = mutation_mask.sum(axis=1, keepdims=True)
    safe_count = np.maximum(count, 1.0)
    mutation_delta = (
        delta * mutation_mask[:, :, None]
    ).sum(axis=1) / safe_count

    wt_rows = count[:, 0] == 0.0
    mutation_delta[wt_rows] = 0.0

    feature = np.concatenate([global_delta, mutation_delta], axis=1)
    if feature.shape[1] != FEATURE_DIM:
        raise RuntimeError(
            f"Unexpected contextual feature dimension: {feature.shape[1]}"
        )
    return feature.astype(np.float32, copy=False)


def encode_context_features(
    sequences: list[str],
    reference: str,
    *,
    batch_size: int,
):
    import torch
    from transformers import AutoModel, AutoTokenizer

    torch.manual_seed(SEED)
    torch.set_num_threads(2)
    try:
        torch.use_deterministic_algorithms(True)
    except Exception:
        pass

    tokenizer = AutoTokenizer.from_pretrained(
        MODEL_NAME,
        revision=MODEL_REVISION,
        trust_remote_code=False,
    )
    encoder = AutoModel.from_pretrained(
        MODEL_NAME,
        revision=MODEL_REVISION,
        trust_remote_code=False,
        use_safetensors=True,
    )
    encoder.eval()
    encoder.to("cpu")

    length = len(reference)

    with torch.inference_mode():
        ref_tokens = tokenizer(
            [reference],
            return_tensors="pt",
            padding=True,
            add_special_tokens=True,
        )
        ref_output = encoder(**ref_tokens).last_hidden_state
        if ref_output.shape[1] < length + 2:
            raise RuntimeError("Unexpected ESM tokenized reference length.")
        ref_residue = (
            ref_output[0, 1 : 1 + length, :]
            .detach()
            .cpu()
            .numpy()
            .astype(np.float32)
        )

    all_features = []
    for start in range(0, len(sequences), batch_size):
        batch = sequences[start : start + batch_size]
        if any(len(sequence) != length for sequence in batch):
            raise RuntimeError("All sequences must match reference length.")

        with torch.inference_mode():
            tokens = tokenizer(
                batch,
                return_tensors="pt",
                padding=True,
                add_special_tokens=True,
            )
            output = encoder(**tokens).last_hidden_state
            residue = (
                output[:, 1 : 1 + length, :]
                .detach()
                .cpu()
                .numpy()
                .astype(np.float32)
            )

        mask = np.asarray(
            [
                [1.0 if aa != ref_aa else 0.0 for aa, ref_aa in zip(sequence, reference)]
                for sequence in batch
            ],
            dtype=np.float32,
        )
        all_features.append(
            contextual_features_from_hidden(
                residue,
                ref_residue,
                mask,
            )
        )

    result = np.concatenate(all_features, axis=0)
    if result.shape != (len(sequences), FEATURE_DIM):
        raise RuntimeError(
            f"Unexpected encoded feature shape: {result.shape}"
        )
    return result


def build_clean_additive_baseline(
    mutation_sets,
    *,
    wt_target,
    clean_main,
    b2_prediction,
):
    baseline = []
    supported = []
    for mutations, fallback in zip(mutation_sets, b2_prediction):
        is_supported = all(
            mutation in clean_main
            for mutation in mutations
        )
        supported.append(is_supported)
        if is_supported:
            baseline.append(
                float(
                    wt_target
                    + sum(clean_main[mutation] for mutation in mutations)
                )
            )
        else:
            baseline.append(float(fallback))
    return np.asarray(baseline, dtype=float), np.asarray(supported, dtype=bool)


def run(source_gz: Path, output_dir: Path, batch_size: int):
    source_sha, fit_frame, validation_frame, test_frame = load_ired(source_gz)
    reference = derive_reference(fit_frame["sequence"].tolist())

    fit_mutations = [
        mutation_set(sequence, reference)
        for sequence in fit_frame["sequence"]
    ]
    test_mutations = [
        mutation_set(sequence, reference)
        for sequence in test_frame["sequence"]
    ]
    fit_target = fit_frame["target"].to_numpy(dtype=float)
    test_target = test_frame["target"].to_numpy(dtype=float)
    test_ids = test_frame["sequence"].astype(str).tolist()

    v83_model = NabuV83Model().fit(
        mutation_sets=fit_mutations,
        labels=fit_target,
        candidate_ids=fit_frame["sequence"].astype(str).tolist(),
    )
    phase2_prediction, strict_scoreable = phase2_abstention_predictions(
        v83_model,
        test_mutations,
        test_ids,
    )
    b2_prediction, b2_all_main = b2_predictions(
        v83_model,
        test_mutations,
    )

    old_factor = fit_pair_factorization(
        v83_model,
        fit_mutations,
        fit_target,
        permute_targets=False,
    )
    v91_prediction = np.asarray(
        [
            score_v91(mutations, v83_model, old_factor)[0]
            for mutations in test_mutations
        ],
        dtype=float,
    )

    clean_metadata, clean = build_clean_pair_table(source_gz)
    wt_target, singleton_target, clean_main = build_singleton_map(
        fit_frame,
        reference,
    )

    old_calibration, _ = fit_global_oof_calibration(clean)
    old_relation = fit_node_additive(
        list(zip(clean["left"], clean["right"])),
        clean["clean_pair_residual"].to_numpy(dtype=float),
        permute=False,
        fold=0,
    )
    old_exact = exact_pair_map(clean)
    old_rows = [
        score_v9_clean_candidate(
            mutations,
            wt_target=wt_target,
            clean_main=clean_main,
            exact_pairs=old_exact,
            relation_model=old_relation,
            calibration=old_calibration,
            b2_fallback=b2_value,
        )
        for mutations, b2_value in zip(test_mutations, b2_prediction)
    ]
    old_clean_prediction = np.asarray(
        [value for value, _ in old_rows],
        dtype=float,
    )

    residual_train = build_residual_training_table(
        fit_frame,
        reference,
        clean,
    )
    train_sequences = residual_train["sequence"].astype(str).tolist()
    train_residual = residual_train[
        "interaction_residual"
    ].to_numpy(dtype=float)

    train_features = encode_context_features(
        train_sequences,
        reference,
        batch_size=batch_size,
    )
    test_features = encode_context_features(
        test_ids,
        reference,
        batch_size=batch_size,
    )

    context_model, train_oof_prediction = fit_context_readout(
        train_features,
        train_residual,
    )
    selected_alpha = float(
        context_model.named_steps["ridge"].alpha_
    )

    contextual_residual = np.asarray(
        context_model.predict(test_features),
        dtype=float,
    )

    clean_baseline, complete_clean_main = build_clean_additive_baseline(
        test_mutations,
        wt_target=wt_target,
        clean_main=clean_main,
        b2_prediction=b2_prediction,
    )
    contextual_prediction = np.where(
        complete_clean_main,
        clean_baseline + contextual_residual,
        b2_prediction,
    ).astype(float)

    replay_residual = np.asarray(
        context_model.predict(test_features),
        dtype=float,
    )
    replay_prediction = np.where(
        complete_clean_main,
        clean_baseline + replay_residual,
        b2_prediction,
    ).astype(float)

    arms = {
        "PHASE2_V8_3_ABSTENTION": phase2_prediction,
        "B2_MAIN_ONLY": b2_prediction,
        "V9_1_DEGREE_NORMALIZED_TRANSFER": v91_prediction,
        "V9_CLEAN_CALIBRATED_PAIR_EXPANSION": old_clean_prediction,
        "V9_CONTEXTUAL_ESM_RESIDUAL": contextual_prediction,
    }
    if not all(np.isfinite(values).all() for values in arms.values()):
        raise RuntimeError("Non-finite prediction detected.")

    metrics = {
        name: evaluate(test_target, values)
        for name, values in arms.items()
    }
    by_mutation_count = {
        name: evaluate_by_mutation_count(
            test_target,
            values,
            test_mutations,
        )
        for name, values in arms.items()
    }

    clean_region_metrics = {
        name: evaluate(
            test_target[complete_clean_main],
            values[complete_clean_main],
        )
        for name, values in {
            "B2_MAIN_ONLY": b2_prediction,
            "V9_CLEAN_CALIBRATED_PAIR_EXPANSION": old_clean_prediction,
            "V9_CONTEXTUAL_ESM_RESIDUAL": contextual_prediction,
        }.items()
    }

    train_oof_all = {
        "count": int(len(train_residual)),
        "spearman": finite_spearman(
            train_residual,
            train_oof_prediction,
        ),
        "rmse": float(
            np.sqrt(
                np.mean(
                    (train_oof_prediction - train_residual) ** 2
                )
            )
        ),
    }
    clean_double_train_mask = (
        residual_train["evidence_class"].to_numpy()
        == "CLEAN_DOUBLE"
    )
    train_oof_clean_double = {
        "count": int(clean_double_train_mask.sum()),
        "spearman": finite_spearman(
            train_residual[clean_double_train_mask],
            train_oof_prediction[clean_double_train_mask],
        ),
        "rmse": float(
            np.sqrt(
                np.mean(
                    (
                        train_oof_prediction[clean_double_train_mask]
                        - train_residual[clean_double_train_mask]
                    )
                    ** 2
                )
            )
        ),
    }

    contextual = metrics["V9_CONTEXTUAL_ESM_RESIDUAL"]
    b2 = metrics["B2_MAIN_ONLY"]
    old_clean = metrics["V9_CLEAN_CALIBRATED_PAIR_EXPANSION"]
    count_metrics = by_mutation_count["V9_CONTEXTUAL_ESM_RESIDUAL"]

    checks = {
        "full_spearman_beats_b2": bool(
            contextual["spearman"] is not None
            and b2["spearman"] is not None
            and contextual["spearman"] > b2["spearman"]
        ),
        "clean_region_spearman_beats_b2": bool(
            clean_region_metrics["V9_CONTEXTUAL_ESM_RESIDUAL"][
                "spearman"
            ]
            > clean_region_metrics["B2_MAIN_ONLY"]["spearman"]
        ),
        "top1_hits_at_least_b2": bool(
            contextual["top1_percent_hits"]
            >= b2["top1_percent_hits"]
        ),
        "regret_not_worse_than_b2": bool(
            contextual["normalized_regret_top1pct"]
            <= b2["normalized_regret_top1pct"]
        ),
        "full_spearman_beats_clean_pair_expansion": bool(
            contextual["spearman"] > old_clean["spearman"]
        ),
        "positive_spearman_counts_3_4_5": bool(
            all(
                str(count) in count_metrics
                and count_metrics[str(count)]["spearman"] is not None
                and count_metrics[str(count)]["spearman"] > 0.0
                for count in (3, 4, 5)
            )
        ),
        "residual_training_oof_spearman_positive": bool(
            train_oof_all["spearman"] is not None
            and train_oof_all["spearman"] > 0.0
        ),
        "deterministic_rank_replay": bool(
            rank_hash(test_ids, contextual_prediction)
            == rank_hash(test_ids, replay_prediction)
        ),
    }

    predictions = pd.DataFrame(
        {
            "candidate_id": test_ids,
            "mutation_count": [len(x) for x in test_mutations],
            "target": test_target,
            "strict_v83_scoreable": strict_scoreable,
            "b2_all_main_supported": b2_all_main,
            "complete_clean_main": complete_clean_main,
            "clean_additive_baseline_or_b2_fallback": clean_baseline,
            "contextual_residual": contextual_residual,
            **arms,
        }
    )

    result = {
        "version": VERSION,
        "status": "POST_PHASE2_DEVELOPMENT_DIAGNOSTIC",
        "phase3_opened": False,
        "scientific_claim": False,
        "source_sha256": source_sha,
        "frozen_v83_core": (
            "c8afdcd7698231d95a39ef13a3fe22b6e3f507c5"
        ),
        "encoder": {
            "model": MODEL_NAME,
            "revision": MODEL_REVISION,
            "frozen": True,
            "feature_definition": (
                "concat(global mean final-hidden delta vs WT, "
                "mean mutation-site final-hidden delta vs WT)"
            ),
            "feature_dim": FEATURE_DIM,
            "train_feature_sha256_rounded_8dp": sha256_array(
                train_features,
                decimals=8,
            ),
            "test_feature_sha256_rounded_8dp": sha256_array(
                test_features,
                decimals=8,
            ),
        },
        "split_counts": {
            "fit": int(len(fit_frame)),
            "validation_unused_for_fit": int(len(validation_frame)),
            "test_development": int(len(test_frame)),
        },
        "clean_evidence": {
            "wt_target": float(wt_target),
            "singleton_count": int(len(singleton_target)),
            "clean_pair_count": int(len(clean)),
            "excluded_double_count": int(
                residual_train.attrs["excluded_double_count"]
            ),
            "residual_training_count": int(len(residual_train)),
            "residual_training_class_counts": {
                str(key): int(value)
                for key, value in residual_train[
                    "evidence_class"
                ].value_counts().to_dict().items()
            },
        },
        "context_readout": {
            "ridge_alpha_grid": RIDGE_ALPHAS.tolist(),
            "selected_alpha": selected_alpha,
            "training_oof_all": train_oof_all,
            "training_oof_clean_doubles": train_oof_clean_double,
        },
        "coverage": {
            "complete_clean_main": int(
                complete_clean_main.sum()
            ),
            "b2_fallback_rows": int(
                (~complete_clean_main).sum()
            ),
        },
        "metrics": metrics,
        "complete_clean_main_region_metrics": clean_region_metrics,
        "by_mutation_count": by_mutation_count,
        "prediction_scale": {
            name: summarize_scale(values)
            for name, values in arms.items()
        },
        "contextual_residual_scale": summarize_scale(
            contextual_residual[complete_clean_main]
        ),
        "prediction_rank_sha256": {
            name: rank_hash(test_ids, values)
            for name, values in arms.items()
        },
        "replay_rank_sha256": rank_hash(
            test_ids,
            replay_prediction,
        ),
        "decision_checks": checks,
        "development_decision": (
            "V9_CONTEXTUAL_ESM_RESIDUAL_RETAIN"
            if all(checks.values())
            else "V9_CONTEXTUAL_ESM_RESIDUAL_REJECT"
        ),
        "interpretation_boundary": (
            "IRED test labels were already revealed in Phase 2D. "
            "This is development evidence only; a retained architecture "
            "requires a new untouched external benchmark."
        ),
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    predictions.to_csv(
        output_dir / "V9_CONTEXTUAL_PREDICTIONS.csv",
        index=False,
    )
    residual_train.assign(
        oof_context_prediction=train_oof_prediction
    ).to_csv(
        output_dir / "V9_CONTEXTUAL_TRAIN_OOF.csv",
        index=False,
    )
    (output_dir / "V9_CONTEXTUAL_RESULTS.json").write_text(
        json.dumps(result, indent=2),
        encoding="utf-8",
    )

    print(json.dumps(result, indent=2))
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source_gz")
    parser.add_argument("--out", required=True)
    parser.add_argument("--batch-size", type=int, default=16)
    args = parser.parse_args()
    run(
        Path(args.source_gz),
        Path(args.out),
        batch_size=args.batch_size,
    )


if __name__ == "__main__":
    main()
