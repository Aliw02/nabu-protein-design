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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("test_masked_csv_gz")
    parser.add_argument("--state", required=True)
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--shard-count", type=int, required=True)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    frame = pd.read_csv(args.test_masked_csv_gz, compression="gzip")
    positions = np.arange(len(frame), dtype=np.int64)
    mask = (positions % args.shard_count) == args.shard_index
    shard = frame.loc[mask].copy().reset_index(drop=True)

    state = joblib.load(args.state)
    reference = state["reference"]
    base = state["base"]
    context_model = state["context_model"]
    shuffled_model = state["shuffled_model"]

    sequences = shard["sequence"].astype(str).tolist()
    mutation_sets = [mutation_set(sequence, reference) for sequence in sequences]

    b2_values = []
    supported = []
    for mutations in mutation_sets:
        value, ok = b2_score(mutations, base)
        b2_values.append(value)
        supported.append(ok)
    b2_values = np.asarray(b2_values, dtype=float)

    features = encode_context_features(
        sequences,
        reference,
        batch_size=args.batch_size,
    )
    context_residual = np.asarray(
        context_model.predict(features),
        dtype=float,
    )
    shuffled_residual = np.asarray(
        shuffled_model.predict(features),
        dtype=float,
    )

    corrected = b2_values + context_residual
    shuffled_corrected = b2_values + shuffled_residual

    for values in (b2_values, corrected, shuffled_corrected):
        if not np.isfinite(values).all():
            raise RuntimeError("Non-finite prediction.")

    output = pd.DataFrame(
        {
            "row_id": shard["row_id"].astype(str),
            "mutation_count": [len(x) for x in mutation_sets],
            "b2_supported": np.asarray(supported, dtype=bool),
            "B2_MAIN_ONLY": b2_values,
            "B2_PLUS_CONTEXT_RESIDUAL": corrected,
            "B2_PLUS_SHUFFLED_CONTEXT": shuffled_corrected,
            "context_residual": context_residual,
            "shuffled_context_residual": shuffled_residual,
        }
    )

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    output.to_csv(
        out / f"pred_{args.shard_index:02d}.csv.gz",
        index=False,
        compression="gzip",
    )

    diag = {
        "shard_index": int(args.shard_index),
        "row_count": int(len(output)),
        "b2_supported_count": int(np.sum(supported)),
    }
    (out / f"diag_{args.shard_index:02d}.json").write_text(
        json.dumps(diag, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(diag, indent=2))


if __name__ == "__main__":
    main()
