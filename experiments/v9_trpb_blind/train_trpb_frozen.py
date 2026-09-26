from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DEV = HERE.parent / "v9_pair_transfer_dev"
sys.path.insert(0, str(DEV))

from nabu_protein.v83 import NabuV83Model
from run_v9_pair_transfer_ired import derive_reference, mutation_set
from run_v9_contextual import encode_context_features, fit_context_readout


def dedup_encode(sequences, reference, batch_size):
    sequences = np.asarray(list(sequences), dtype=str)
    unique, inverse = np.unique(sequences, return_inverse=True)
    unique_features = encode_context_features(
        unique.tolist(),
        reference,
        batch_size=batch_size,
    )
    return unique_features[inverse], int(len(unique))


def finite_spearman(target, prediction):
    from scipy.stats import spearmanr

    value = float(spearmanr(target, prediction).statistic)
    return value if np.isfinite(value) else None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fit_csv_gz")
    parser.add_argument("--source-manifest", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--batch-size", type=int, default=32)
    args = parser.parse_args()

    fit = pd.read_csv(args.fit_csv_gz, compression="gzip")
    source_manifest = json.loads(
        Path(args.source_manifest).read_text(encoding="utf-8")
    )
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    reference = derive_reference(fit["sequence"].astype(str).tolist())
    mutations = [
        mutation_set(sequence, reference)
        for sequence in fit["sequence"].astype(str)
    ]
    target = fit["target"].to_numpy(dtype=float)
    row_ids = fit["row_id"].astype(str).tolist()

    model = NabuV83Model().fit(
        mutation_sets=mutations,
        labels=target,
        candidate_ids=row_ids,
    )
    base = model.model["base"]

    wt_targets = []
    singleton_targets = defaultdict(list)
    for muts, value in zip(mutations, target):
        if len(muts) == 0:
            wt_targets.append(float(value))
        elif len(muts) == 1:
            singleton_targets[muts[0]].append(float(value))

    if not wt_targets:
        raise RuntimeError("No WT evidence in TrpB fit/train.")

    wt_target = float(np.mean(wt_targets))
    singleton_mean = {
        mutation: float(np.mean(values))
        for mutation, values in singleton_targets.items()
    }
    clean_main = {
        mutation: float(value - wt_target)
        for mutation, value in singleton_mean.items()
    }

    residual_sequences = []
    residual_targets = []
    residual_classes = []
    excluded_double = 0

    for sequence, muts, value in zip(
        fit["sequence"].astype(str),
        mutations,
        target,
    ):
        if len(muts) == 0:
            residual_sequences.append(sequence)
            residual_targets.append(0.0)
            residual_classes.append("WT")
        elif len(muts) == 1:
            residual_sequences.append(sequence)
            residual_targets.append(0.0)
            residual_classes.append("SINGLE")
        elif len(muts) == 2:
            left, right = muts
            if left not in singleton_mean or right not in singleton_mean:
                excluded_double += 1
                continue
            residual = (
                float(value)
                - singleton_mean[left]
                - singleton_mean[right]
                + wt_target
            )
            residual_sequences.append(sequence)
            residual_targets.append(float(residual))
            residual_classes.append("CLEAN_DOUBLE")
        else:
            raise RuntimeError("Fit/train mutation order exceeded 2.")

    residual_targets = np.asarray(residual_targets, dtype=float)
    features, unique_encoded = dedup_encode(
        residual_sequences,
        reference,
        args.batch_size,
    )
    context_model, oof_prediction = fit_context_readout(
        features,
        residual_targets,
    )

    clean_double_mask = np.asarray(
        [value == "CLEAN_DOUBLE" for value in residual_classes],
        dtype=bool,
    )
    if int(clean_double_mask.sum()) == 0:
        raise RuntimeError("No clean double residual evidence.")

    all_rmse = float(
        np.sqrt(np.mean((oof_prediction - residual_targets) ** 2))
    )
    clean_rmse = float(
        np.sqrt(
            np.mean(
                (
                    oof_prediction[clean_double_mask]
                    - residual_targets[clean_double_mask]
                )
                ** 2
            )
        )
    )

    state = {
        "version": "NABU_V9_TRPB_BLIND_FROZEN_STATE_V1",
        "reference": reference,
        "source_sha256": source_manifest["source_sha256"],
        "base": base,
        "wt_target": wt_target,
        "singleton_mean": singleton_mean,
        "clean_main": clean_main,
        "context_model": context_model,
        "validation_targets_used_for_fit": False,
        "architecture": "V9_CONSENSUS_B2_CONTEXT",
    }
    state_path = out / "TRPB_FROZEN_STATE.joblib"
    joblib.dump(state, state_path)

    state_sha = hashlib.sha256(state_path.read_bytes()).hexdigest()
    diagnostics = {
        "fit_count": int(len(fit)),
        "wt_row_count": int(len(wt_targets)),
        "singleton_identity_count": int(len(singleton_mean)),
        "context_residual_training_count": int(len(residual_targets)),
        "context_residual_unique_sequence_count": unique_encoded,
        "excluded_double_count": int(excluded_double),
        "clean_double_count": int(clean_double_mask.sum()),
        "selected_ridge_alpha": float(
            context_model.named_steps["ridge"].alpha_
        ),
        "oof_all_spearman": finite_spearman(
            residual_targets,
            oof_prediction,
        ),
        "oof_all_rmse": all_rmse,
        "oof_clean_double_spearman": finite_spearman(
            residual_targets[clean_double_mask],
            oof_prediction[clean_double_mask],
        ),
        "oof_clean_double_rmse": clean_rmse,
        "state_sha256": state_sha,
        "validation_targets_used_for_fit": False,
    }
    (out / "TRPB_TRAIN_DIAGNOSTICS.json").write_text(
        json.dumps(diagnostics, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(diagnostics, indent=2))


if __name__ == "__main__":
    main()
