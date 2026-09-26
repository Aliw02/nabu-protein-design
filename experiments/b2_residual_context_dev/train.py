from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
DEV = HERE.parent / "v9_pair_transfer_dev"
sys.path.insert(0, str(DEV))

from nabu_protein.v83 import NabuV83Model
from run_v9_pair_transfer_ired import derive_reference, mutation_set
from run_v9_contextual import encode_context_features, fit_context_readout

SEED = 161


def finite_spearman(a, b):
    value = float(spearmanr(a, b).statistic)
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

    sequences = fit["sequence"].astype(str).str.strip().str.upper().tolist()
    target = fit["target"].to_numpy(dtype=float)
    ids = fit["row_id"].astype(str).tolist()

    reference = derive_reference(sequences)
    mutation_sets = [mutation_set(sequence, reference) for sequence in sequences]

    model = NabuV83Model().fit(
        mutation_sets=mutation_sets,
        labels=target,
        candidate_ids=ids,
    )

    oof_b2 = np.asarray(model.model["oof_b2"], dtype=float)
    if len(oof_b2) != len(target) or not np.isfinite(oof_b2).all():
        raise RuntimeError("Frozen V8.3 OOF B2 predictions unavailable or non-finite.")

    residual_target = target - oof_b2

    features = encode_context_features(
        sequences,
        reference,
        batch_size=args.batch_size,
    )
    context_model, context_oof = fit_context_readout(
        features,
        residual_target,
    )

    rng = np.random.default_rng(SEED)
    shuffled_target = rng.permutation(residual_target)
    shuffled_model, shuffled_oof = fit_context_readout(
        features,
        shuffled_target,
    )

    corrected_oof = oof_b2 + context_oof
    shuffled_corrected_oof = oof_b2 + shuffled_oof

    state = {
        "version": "NABU_B2_RESIDUAL_CONTEXT_DEV_STATE_V1",
        "reference": reference,
        "source_sha256": source_manifest["source_sha256"],
        "base": model.model["base"],
        "context_model": context_model,
        "shuffled_model": shuffled_model,
        "architecture": "B2_PLUS_CONTEXT_RESIDUAL",
        "validation_targets_used_for_fit": False,
    }
    state_path = out / "B2_RESIDUAL_CONTEXT_STATE.joblib"
    joblib.dump(state, state_path)

    state_sha = hashlib.sha256(state_path.read_bytes()).hexdigest()
    diagnostics = {
        "fit_count": int(len(fit)),
        "reference_length": int(len(reference)),
        "mutation_order_distribution": {
            str(order): int(sum(len(m) == order for m in mutation_sets))
            for order in sorted(set(len(m) for m in mutation_sets))
        },
        "selected_ridge_alpha": float(
            context_model.named_steps["ridge"].alpha_
        ),
        "selected_shuffled_ridge_alpha": float(
            shuffled_model.named_steps["ridge"].alpha_
        ),
        "residual_target_mean": float(np.mean(residual_target)),
        "residual_target_std": float(np.std(residual_target)),
        "b2_oof_spearman": finite_spearman(target, oof_b2),
        "corrected_oof_spearman": finite_spearman(target, corrected_oof),
        "shuffled_corrected_oof_spearman": finite_spearman(
            target, shuffled_corrected_oof
        ),
        "context_residual_oof_spearman": finite_spearman(
            residual_target, context_oof
        ),
        "shuffled_residual_oof_spearman": finite_spearman(
            residual_target, shuffled_oof
        ),
        "state_sha256": state_sha,
        "validation_targets_used_for_fit": False,
    }
    (out / "TRAIN_DIAGNOSTICS.json").write_text(
        json.dumps(diagnostics, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(diagnostics, indent=2))


if __name__ == "__main__":
    main()
