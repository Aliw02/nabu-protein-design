from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr
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
)
from run_v9_contextual import (
    MODEL_NAME,
    MODEL_REVISION,
    RIDGE_ALPHAS,
    SEED,
    build_clean_additive_baseline,
    build_residual_training_table,
)
from run_clean_pair_poc import build_clean_pair_table
from run_v9_candidate import build_singleton_map


VERSION = "NABU_V9_SIMPLE_ARCHITECTURE_ABLATION_V1"
BASE_FEATURE_DIM = 640
CONTACT_FEATURE_DIM = 320
CONTACT_AWARE_FEATURE_DIM = 960


def finite_spearman(target, prediction):
    value = float(spearmanr(target, prediction).statistic)
    return value if np.isfinite(value) else None


def sha256_array(values: np.ndarray, decimals: int = 8) -> str:
    array = np.asarray(values, dtype=np.float64)
    rounded = np.round(array, decimals=decimals)
    return hashlib.sha256(rounded.tobytes(order="C")).hexdigest()


def rank_hash(candidate_ids, prediction) -> str:
    frame = pd.DataFrame(
        {
            "candidate_id": list(candidate_ids),
            "prediction": np.asarray(prediction, dtype=float),
        }
    )
    order = (
        frame.sort_values(
            ["prediction", "candidate_id"],
            ascending=[False, True],
            kind="stable",
        )["candidate_id"]
        .astype(str)
        .tolist()
    )
    return hashlib.sha256("\n".join(order).encode("utf-8")).hexdigest()


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


def fit_readout(features: np.ndarray, target: np.ndarray):
    features = np.asarray(features, dtype=float)
    target = np.asarray(target, dtype=float)
    if features.ndim != 2 or len(features) != len(target):
        raise ValueError("Feature/target shape mismatch.")
    if not np.isfinite(features).all() or not np.isfinite(target).all():
        raise ValueError("Readout requires finite inputs.")

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


def base_context_features(
    residue_hidden: np.ndarray,
    reference_hidden: np.ndarray,
    mutation_mask: np.ndarray,
) -> np.ndarray:
    delta = residue_hidden - reference_hidden[None, :, :]
    global_delta = delta.mean(axis=1)

    mutation_count = mutation_mask.sum(axis=1, keepdims=True)
    safe_count = np.maximum(mutation_count, 1.0)
    mutation_delta = (
        delta * mutation_mask[:, :, None]
    ).sum(axis=1) / safe_count
    mutation_delta[mutation_count[:, 0] == 0.0] = 0.0

    feature = np.concatenate([global_delta, mutation_delta], axis=1)
    if feature.shape[1] != BASE_FEATURE_DIM:
        raise RuntimeError(f"Unexpected base feature shape: {feature.shape}")
    return feature.astype(np.float32, copy=False)


def contact_neighbor_features(
    residue_hidden: np.ndarray,
    reference_hidden: np.ndarray,
    mutation_mask: np.ndarray,
    contacts: np.ndarray,
) -> np.ndarray:
    residue_hidden = np.asarray(residue_hidden, dtype=np.float32)
    reference_hidden = np.asarray(reference_hidden, dtype=np.float32)
    mutation_mask = np.asarray(mutation_mask, dtype=np.float32)
    contacts = np.asarray(contacts, dtype=np.float32)

    batch, length, hidden = residue_hidden.shape
    if reference_hidden.shape != (length, hidden):
        raise ValueError("Reference hidden shape mismatch.")
    if mutation_mask.shape != (batch, length):
        raise ValueError("Mutation mask shape mismatch.")
    if contacts.shape != (batch, length, length):
        raise ValueError(
            f"Contact shape mismatch: {contacts.shape} vs {(batch, length, length)}"
        )

    delta = residue_hidden - reference_hidden[None, :, :]

    weights = contacts.copy()
    diagonal = np.arange(length)
    weights[:, diagonal, diagonal] = 0.0

    denominator = weights.sum(axis=2, keepdims=True)
    denominator = np.maximum(denominator, 1e-6)
    neighborhood = np.einsum(
        "bij,bjh->bih",
        weights,
        delta,
        optimize=True,
    ) / denominator

    mutation_count = mutation_mask.sum(axis=1, keepdims=True)
    safe_count = np.maximum(mutation_count, 1.0)
    contact_context = (
        neighborhood * mutation_mask[:, :, None]
    ).sum(axis=1) / safe_count

    contact_context[mutation_count[:, 0] == 0.0] = 0.0
    if contact_context.shape[1] != CONTACT_FEATURE_DIM:
        raise RuntimeError(
            f"Unexpected contact feature shape: {contact_context.shape}"
        )
    return contact_context.astype(np.float32, copy=False)


def encode_base_and_contact_features(
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
        ref_hidden = encoder(**ref_tokens).last_hidden_state
        ref_residue = (
            ref_hidden[0, 1 : 1 + length, :]
            .detach()
            .cpu()
            .numpy()
            .astype(np.float32)
        )

    base_chunks = []
    contact_chunks = []

    for start in range(0, len(sequences), batch_size):
        batch = sequences[start : start + batch_size]
        if any(len(sequence) != length for sequence in batch):
            raise RuntimeError("Sequence/reference length mismatch.")

        with torch.inference_mode():
            tokens = tokenizer(
                batch,
                return_tensors="pt",
                padding=True,
                add_special_tokens=True,
            )
            output = encoder(**tokens).last_hidden_state
            contacts = encoder.predict_contacts(
                tokens["input_ids"],
                tokens["attention_mask"],
            )

        residue = (
            output[:, 1 : 1 + length, :]
            .detach()
            .cpu()
            .numpy()
            .astype(np.float32)
        )
        contact_array = (
            contacts[:, :length, :length]
            .detach()
            .cpu()
            .numpy()
            .astype(np.float32)
        )
        mask = np.asarray(
            [
                [
                    1.0 if aa != ref_aa else 0.0
                    for aa, ref_aa in zip(sequence, reference)
                ]
                for sequence in batch
            ],
            dtype=np.float32,
        )

        base_chunks.append(
            base_context_features(
                residue,
                ref_residue,
                mask,
            )
        )
        contact_chunks.append(
            contact_neighbor_features(
                residue,
                ref_residue,
                mask,
                contact_array,
            )
        )

    base = np.concatenate(base_chunks, axis=0)
    contact = np.concatenate(contact_chunks, axis=0)
    combined = np.concatenate([base, contact], axis=1)

    if base.shape != (len(sequences), BASE_FEATURE_DIM):
        raise RuntimeError(f"Unexpected base matrix shape: {base.shape}")
    if contact.shape != (len(sequences), CONTACT_FEATURE_DIM):
        raise RuntimeError(f"Unexpected contact matrix shape: {contact.shape}")
    if combined.shape != (len(sequences), CONTACT_AWARE_FEATURE_DIM):
        raise RuntimeError(f"Unexpected combined shape: {combined.shape}")

    return base, contact, combined


def pareto_consensus_score(
    b2_score: np.ndarray,
    context_score: np.ndarray,
    candidate_ids: list[str],
) -> tuple[np.ndarray, np.ndarray]:
    b2_score = np.asarray(b2_score, dtype=float)
    context_score = np.asarray(context_score, dtype=float)
    n = len(b2_score)

    if len(context_score) != n or len(candidate_ids) != n:
        raise ValueError("Consensus input length mismatch.")

    b2_rank = rankdata(-b2_score, method="average") / float(n)
    context_rank = rankdata(-context_score, method="average") / float(n)

    dominates = (
        (b2_rank[:, None] <= b2_rank[None, :])
        & (context_rank[:, None] <= context_rank[None, :])
        & (
            (b2_rank[:, None] < b2_rank[None, :])
            | (context_rank[:, None] < context_rank[None, :])
        )
    )
    dominator_count = dominates.sum(axis=0).astype(np.int64)

    layer = np.full(n, -1, dtype=np.int64)
    remaining = np.ones(n, dtype=bool)
    current_layer = 0

    while remaining.any():
        front = np.where(remaining & (dominator_count == 0))[0]
        if len(front) == 0:
            raise RuntimeError("Pareto sorting reached an empty front.")
        layer[front] = current_layer
        remaining[front] = False

        decrement = dominates[front].sum(axis=0).astype(np.int64)
        dominator_count[remaining] -= decrement[remaining]
        current_layer += 1

    mean_rank = 0.5 * (b2_rank + context_rank)
    ids = np.asarray(candidate_ids, dtype=str)
    order = np.lexsort((ids, mean_rank, layer))

    score = np.empty(n, dtype=float)
    score[order] = np.arange(n, 0, -1, dtype=float)
    return score, layer


def train_oof_metrics(target, prediction):
    target = np.asarray(target, dtype=float)
    prediction = np.asarray(prediction, dtype=float)
    return {
        "count": int(len(target)),
        "spearman": finite_spearman(target, prediction),
        "rmse": float(np.sqrt(np.mean((prediction - target) ** 2))),
    }


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
    b2_prediction, _ = b2_predictions(v83_model, test_mutations)

    clean_metadata, clean = build_clean_pair_table(source_gz)
    wt_target, singleton_target, clean_main = build_singleton_map(
        fit_frame,
        reference,
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

    train_base, train_contact, train_combined = encode_base_and_contact_features(
        train_sequences,
        reference,
        batch_size=batch_size,
    )
    test_base, test_contact, test_combined = encode_base_and_contact_features(
        test_ids,
        reference,
        batch_size=batch_size,
    )

    context_model, context_oof = fit_readout(
        train_base,
        train_residual,
    )
    contact_model, contact_oof = fit_readout(
        train_combined,
        train_residual,
    )

    context_residual = np.asarray(
        context_model.predict(test_base),
        dtype=float,
    )
    contact_residual = np.asarray(
        contact_model.predict(test_combined),
        dtype=float,
    )

    clean_baseline, complete_clean_main = build_clean_additive_baseline(
        test_mutations,
        wt_target=wt_target,
        clean_main=clean_main,
        b2_prediction=b2_prediction,
    )

    context_prediction = np.where(
        complete_clean_main,
        clean_baseline + context_residual,
        b2_prediction,
    ).astype(float)

    contact_prediction = np.where(
        complete_clean_main,
        clean_baseline + contact_residual,
        b2_prediction,
    ).astype(float)

    consensus_prediction, pareto_layer = pareto_consensus_score(
        b2_prediction,
        context_prediction,
        test_ids,
    )

    consensus_replay, replay_layer = pareto_consensus_score(
        b2_prediction,
        context_prediction,
        test_ids,
    )
    contact_replay = np.where(
        complete_clean_main,
        clean_baseline
        + np.asarray(contact_model.predict(test_combined), dtype=float),
        b2_prediction,
    ).astype(float)

    arms = {
        "B2_MAIN_ONLY": b2_prediction,
        "V9_CONTEXTUAL_ESM_RESIDUAL": context_prediction,
        "V9_CONSENSUS_B2_CONTEXT": consensus_prediction,
        "V9_CONTACT_AWARE_CONTEXT": contact_prediction,
    }

    if not all(np.isfinite(values).all() for values in arms.values()):
        raise RuntimeError("Non-finite ablation prediction.")

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
        for name, values in arms.items()
    }

    b2 = metrics["B2_MAIN_ONLY"]

    def common_checks(arm_name: str, replay_prediction: np.ndarray):
        arm = metrics[arm_name]
        clean_arm = clean_region_metrics[arm_name]
        clean_b2 = clean_region_metrics["B2_MAIN_ONLY"]
        order_metrics = by_mutation_count[arm_name]
        return {
            "full_spearman_at_least_b2": bool(
                arm["spearman"] is not None
                and b2["spearman"] is not None
                and arm["spearman"] >= b2["spearman"]
            ),
            "clean_region_spearman_at_least_b2": bool(
                clean_arm["spearman"] is not None
                and clean_b2["spearman"] is not None
                and clean_arm["spearman"] >= clean_b2["spearman"]
            ),
            "top1_hits_at_least_b2": bool(
                arm["top1_percent_hits"] >= b2["top1_percent_hits"]
            ),
            "regret_not_worse_than_b2": bool(
                arm["normalized_regret_top1pct"]
                <= b2["normalized_regret_top1pct"]
            ),
            "positive_spearman_counts_3_4_5": bool(
                all(
                    str(count) in order_metrics
                    and order_metrics[str(count)]["spearman"] is not None
                    and order_metrics[str(count)]["spearman"] > 0.0
                    for count in (3, 4, 5)
                )
            ),
            "deterministic_replay": bool(
                rank_hash(test_ids, arms[arm_name])
                == rank_hash(test_ids, replay_prediction)
            ),
        }

    consensus_checks = common_checks(
        "V9_CONSENSUS_B2_CONTEXT",
        consensus_replay,
    )
    contact_checks = common_checks(
        "V9_CONTACT_AWARE_CONTEXT",
        contact_replay,
    )
    contact_oof_result = train_oof_metrics(
        train_residual,
        contact_oof,
    )
    contact_checks[
        "residual_training_oof_spearman_positive"
    ] = bool(
        contact_oof_result["spearman"] is not None
        and contact_oof_result["spearman"] > 0.0
    )

    consensus_pass = all(consensus_checks.values())
    contact_pass = all(contact_checks.values())

    if consensus_pass and not contact_pass:
        winner = "V9_CONSENSUS_B2_CONTEXT"
    elif contact_pass and not consensus_pass:
        winner = "V9_CONTACT_AWARE_CONTEXT"
    elif consensus_pass and contact_pass:
        consensus_metric = metrics["V9_CONSENSUS_B2_CONTEXT"]
        contact_metric = metrics["V9_CONTACT_AWARE_CONTEXT"]
        delta = (
            consensus_metric["spearman"]
            - contact_metric["spearman"]
        )
        if abs(delta) > 1e-12:
            winner = (
                "V9_CONSENSUS_B2_CONTEXT"
                if delta > 0
                else "V9_CONTACT_AWARE_CONTEXT"
            )
        else:
            clean_consensus = clean_region_metrics[
                "V9_CONSENSUS_B2_CONTEXT"
            ]["spearman"]
            clean_contact = clean_region_metrics[
                "V9_CONTACT_AWARE_CONTEXT"
            ]["spearman"]
            if abs(clean_consensus - clean_contact) > 1e-12:
                winner = (
                    "V9_CONSENSUS_B2_CONTEXT"
                    if clean_consensus > clean_contact
                    else "V9_CONTACT_AWARE_CONTEXT"
                )
            else:
                winner = (
                    "V9_CONSENSUS_B2_CONTEXT"
                    if consensus_metric["top1_percent_hits"]
                    >= contact_metric["top1_percent_hits"]
                    else "V9_CONTACT_AWARE_CONTEXT"
                )
    else:
        winner = "NO_SIMPLE_SOLUTION_PASS"

    predictions = pd.DataFrame(
        {
            "candidate_id": test_ids,
            "mutation_count": [len(x) for x in test_mutations],
            "target": test_target,
            "complete_clean_main": complete_clean_main,
            "pareto_layer": pareto_layer,
            "contextual_residual": context_residual,
            "contact_residual": contact_residual,
            **arms,
        }
    )

    result = {
        "version": VERSION,
        "status": "POST_PHASE2_DEVELOPMENT_ABLATION",
        "phase3_opened": False,
        "scientific_claim": False,
        "source_sha256": source_sha,
        "split_counts": {
            "fit": int(len(fit_frame)),
            "validation_unused_for_fit": int(len(validation_frame)),
            "test_development": int(len(test_frame)),
        },
        "encoder": {
            "model": MODEL_NAME,
            "revision": MODEL_REVISION,
            "frozen": True,
        },
        "features": {
            "base_context_dim": BASE_FEATURE_DIM,
            "contact_neighbor_dim": CONTACT_FEATURE_DIM,
            "contact_aware_dim": CONTACT_AWARE_FEATURE_DIM,
            "train_base_sha256_rounded_8dp": sha256_array(train_base),
            "train_contact_sha256_rounded_8dp": sha256_array(train_contact),
            "test_base_sha256_rounded_8dp": sha256_array(test_base),
            "test_contact_sha256_rounded_8dp": sha256_array(test_contact),
        },
        "readouts": {
            "context_selected_alpha": float(
                context_model.named_steps["ridge"].alpha_
            ),
            "contact_selected_alpha": float(
                contact_model.named_steps["ridge"].alpha_
            ),
            "context_training_oof": train_oof_metrics(
                train_residual,
                context_oof,
            ),
            "contact_training_oof": contact_oof_result,
        },
        "coverage": {
            "complete_clean_main": int(complete_clean_main.sum()),
            "b2_fallback_rows": int((~complete_clean_main).sum()),
        },
        "metrics": metrics,
        "complete_clean_main_region_metrics": clean_region_metrics,
        "by_mutation_count": by_mutation_count,
        "consensus": {
            "pareto_layer_count": int(pareto_layer.max() + 1),
            "pareto_layer_sha256": hashlib.sha256(
                pareto_layer.astype(np.int64).tobytes()
            ).hexdigest(),
        },
        "rank_sha256": {
            name: rank_hash(test_ids, values)
            for name, values in arms.items()
        },
        "consensus_replay_layer_equal": bool(
            np.array_equal(pareto_layer, replay_layer)
        ),
        "decision_checks": {
            "V9_CONSENSUS_B2_CONTEXT": consensus_checks,
            "V9_CONTACT_AWARE_CONTEXT": contact_checks,
        },
        "arm_pass": {
            "V9_CONSENSUS_B2_CONTEXT": bool(consensus_pass),
            "V9_CONTACT_AWARE_CONTEXT": bool(contact_pass),
        },
        "ablation_winner": winner,
        "interpretation_boundary": (
            "IRED test labels were previously revealed. "
            "This ablation is development evidence only and cannot "
            "repair Phase 2 or open Phase 3."
        ),
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    predictions.to_csv(
        output_dir / "V9_SIMPLE_ABLATION_PREDICTIONS.csv",
        index=False,
    )
    residual_train.assign(
        context_oof=context_oof,
        contact_oof=contact_oof,
    ).to_csv(
        output_dir / "V9_SIMPLE_ABLATION_TRAIN_OOF.csv",
        index=False,
    )
    (output_dir / "V9_SIMPLE_ABLATION_RESULTS.json").write_text(
        json.dumps(result, indent=2),
        encoding="utf-8",
    )

    print(json.dumps(result, indent=2))
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source_gz")
    parser.add_argument("--out", required=True)
    parser.add_argument("--batch-size", type=int, default=8)
    args = parser.parse_args()
    run(
        Path(args.source_gz),
        Path(args.out),
        batch_size=args.batch_size,
    )


if __name__ == "__main__":
    main()
