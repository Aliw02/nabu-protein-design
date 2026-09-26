from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from scalable_consensus import (
    pareto_consensus_score_scalable,
    rank_hash,
)


EXPECTED_TEST_COUNT = 217507


def semantic_hash(frame: pd.DataFrame) -> str:
    ordered = frame.sort_values("row_id", kind="stable")
    payload = "\n".join(
        (
            f"{row.row_id},"
            f"{float(row.B2_MAIN_ONLY):.17g},"
            f"{float(row.V9_CONTEXTUAL_ESM_RESIDUAL):.17g},"
            f"{float(row.V9_CONSENSUS_B2_CONTEXT):.17g}"
        )
        for row in ordered.itertuples(index=False)
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("shard_dir")
    parser.add_argument("--source-manifest", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    shard_dir = Path(args.shard_dir)
    source_manifest = json.loads(
        Path(args.source_manifest).read_text(encoding="utf-8")
    )

    paths = sorted(shard_dir.glob("pred_*.csv.gz"))
    if not paths:
        raise RuntimeError("No shard prediction files found.")

    frames = [
        pd.read_csv(path, compression="gzip")
        for path in paths
    ]
    frame = pd.concat(frames, ignore_index=True)

    if len(frame) != EXPECTED_TEST_COUNT:
        raise RuntimeError(
            f"Expected {EXPECTED_TEST_COUNT} test predictions, got {len(frame)}"
        )
    if frame["row_id"].duplicated().any():
        raise RuntimeError("Duplicate row_id in aggregated predictions.")

    frame = frame.sort_values("row_id", kind="stable").reset_index(drop=True)
    candidate_ids = frame["row_id"].astype(str).to_numpy()

    consensus, layer, b2_rank, context_rank = (
        pareto_consensus_score_scalable(
            frame["B2_MAIN_ONLY"].to_numpy(dtype=float),
            frame["V9_CONTEXTUAL_ESM_RESIDUAL"].to_numpy(dtype=float),
            candidate_ids,
        )
    )
    frame["pareto_layer"] = layer
    frame["b2_fractional_rank"] = b2_rank
    frame["context_fractional_rank"] = context_rank
    frame["V9_CONSENSUS_B2_CONTEXT"] = consensus

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    frozen_path = out / "TRPB_FROZEN_PREDICTIONS.csv.gz"
    frame.to_csv(
        frozen_path,
        index=False,
        compression="gzip",
    )

    freeze = {
        "version": "NABU_V9_TRPB_FROZEN_PREDICTIONS_V1",
        "architecture": "V9_CONSENSUS_B2_CONTEXT",
        "source_sha256": source_manifest["source_sha256"],
        "test_count": int(len(frame)),
        "test_targets_visible_during_prediction": False,
        "prediction_semantic_sha256": semantic_hash(frame),
        "b2_rank_sha256": rank_hash(
            candidate_ids,
            frame["B2_MAIN_ONLY"].to_numpy(dtype=float),
        ),
        "context_rank_sha256": rank_hash(
            candidate_ids,
            frame["V9_CONTEXTUAL_ESM_RESIDUAL"].to_numpy(dtype=float),
        ),
        "consensus_rank_sha256": rank_hash(
            candidate_ids,
            consensus,
        ),
        "pareto_layer_sha256": hashlib.sha256(
            layer.astype(np.int64).tobytes()
        ).hexdigest(),
        "pareto_layer_count": int(layer.max() + 1),
        "freeze_completed_before_reveal": True,
    }
    (out / "TRPB_FREEZE_MANIFEST.json").write_text(
        json.dumps(freeze, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(freeze, indent=2))


if __name__ == "__main__":
    main()
