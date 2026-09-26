from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DEV = HERE.parent / "v9_pair_transfer_dev"
sys.path.insert(0, str(DEV))

from nabu_protein.higher_order import additive_score
from run_v9_pair_transfer_ired import mutation_set
from run_v9_contextual import encode_context_features


def b2_score(mutations, base):
    main_memory = base["main"]
    global_mean = float(base["global_mean"])
    if not all(mutation in main_memory for mutation in mutations):
        return global_mean, False
    return (
        float(additive_score(mutations, global_mean, main_memory)),
        True,
    )


def dedup_predict_context(
    sequences,
    reference,
    context_model,
    batch_size,
):
    sequences = np.asarray(list(sequences), dtype=str)
    if len(sequences) == 0:
        return np.asarray([], dtype=float), 0

    unique, inverse = np.unique(sequences, return_inverse=True)
    features = encode_context_features(
        unique.tolist(),
        reference,
        batch_size=batch_size,
    )
    unique_prediction = np.asarray(
        context_model.predict(features),
        dtype=float,
    )
    return unique_prediction[inverse], int(len(unique))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("test_masked_csv_gz")
    parser.add_argument("--state", required=True)
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--shard-count", type=int, required=True)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    if not (0 <= args.shard_index < args.shard_count):
        raise ValueError("Invalid shard index/count.")

    frame = pd.read_csv(args.test_masked_csv_gz, compression="gzip")
    positions = np.arange(len(frame), dtype=np.int64)
    mask = (positions % args.shard_count) == args.shard_index
    shard = frame.loc[mask].copy().reset_index(drop=True)

    state = joblib.load(args.state)
    reference = state["reference"]
    base = state["base"]
    clean_main = state["clean_main"]
    wt_target = float(state["wt_target"])
    context_model = state["context_model"]

    mutation_sets = [
        mutation_set(sequence, reference)
        for sequence in shard["sequence"].astype(str)
    ]

    b2_values = []
    b2_supported = []
    clean_supported = []
    clean_baseline = []

    for mutations in mutation_sets:
        b2_value, supported = b2_score(mutations, base)
        b2_values.append(b2_value)
        b2_supported.append(supported)

        clean_ok = all(mutation in clean_main for mutation in mutations)
        clean_supported.append(clean_ok)
        if clean_ok:
            clean_baseline.append(
                float(
                    wt_target
                    + sum(clean_main[mutation] for mutation in mutations)
                )
            )
        else:
            clean_baseline.append(float(b2_value))

    b2_values = np.asarray(b2_values, dtype=float)
    clean_supported = np.asarray(clean_supported, dtype=bool)
    clean_baseline = np.asarray(clean_baseline, dtype=float)

    context_residual = np.zeros(len(shard), dtype=float)
    supported_indices = np.where(clean_supported)[0]
    if len(supported_indices):
        supported_sequences = (
            shard.iloc[supported_indices]["sequence"].astype(str).tolist()
        )
        residual_values, unique_count = dedup_predict_context(
            supported_sequences,
            reference,
            context_model,
            args.batch_size,
        )
        context_residual[supported_indices] = residual_values
    else:
        unique_count = 0

    context_score = np.where(
        clean_supported,
        clean_baseline + context_residual,
        b2_values,
    ).astype(float)

    if not np.isfinite(b2_values).all():
        raise RuntimeError("Non-finite B2 prediction.")
    if not np.isfinite(context_score).all():
        raise RuntimeError("Non-finite Context prediction.")

    output = pd.DataFrame(
        {
            "row_id": shard["row_id"].astype(str),
            "mutation_count": [len(x) for x in mutation_sets],
            "b2_supported": np.asarray(b2_supported, dtype=bool),
            "complete_clean_main": clean_supported,
            "B2_MAIN_ONLY": b2_values,
            "V9_CONTEXTUAL_ESM_RESIDUAL": context_score,
        }
    )

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"pred_{args.shard_index:02d}.csv.gz"
    output.to_csv(path, index=False, compression="gzip")

    diag = {
        "shard_index": int(args.shard_index),
        "shard_count": int(args.shard_count),
        "row_count": int(len(output)),
        "clean_supported_count": int(clean_supported.sum()),
        "unique_supported_sequences_encoded": int(unique_count),
        "test_targets_visible": False,
    }
    (out / f"diag_{args.shard_index:02d}.json").write_text(
        json.dumps(diag, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(diag, indent=2))


if __name__ == "__main__":
    main()
