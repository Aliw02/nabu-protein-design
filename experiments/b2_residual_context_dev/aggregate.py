from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


def semantic_hash(frame: pd.DataFrame) -> str:
    ordered = frame.sort_values("row_id", kind="stable")
    payload = "\n".join(
        (
            f"{row.row_id},"
            f"{float(row.B2_MAIN_ONLY):.17g},"
            f"{float(row.B2_PLUS_CONTEXT_RESIDUAL):.17g},"
            f"{float(row.B2_PLUS_SHUFFLED_CONTEXT):.17g}"
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

    source = json.loads(
        Path(args.source_manifest).read_text(encoding="utf-8")
    )
    paths = sorted(Path(args.shard_dir).glob("pred_*.csv.gz"))
    if not paths:
        raise RuntimeError("No prediction shards.")

    frame = pd.concat(
        [pd.read_csv(path, compression="gzip") for path in paths],
        ignore_index=True,
    )
    expected = int(source["observed_counts"]["test"])
    if len(frame) != expected:
        raise RuntimeError(f"Expected {expected} rows, got {len(frame)}")
    if frame["row_id"].duplicated().any():
        raise RuntimeError("Duplicate row_id.")

    frame = frame.sort_values("row_id", kind="stable").reset_index(drop=True)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    frozen = out / "FROZEN_PREDICTIONS.csv.gz"
    frame.to_csv(frozen, index=False, compression="gzip")

    manifest = {
        "version": "NABU_B2_RESIDUAL_CONTEXT_DEV_FREEZE_V1",
        "source_sha256": source["source_sha256"],
        "test_count": int(len(frame)),
        "prediction_semantic_sha256": semantic_hash(frame),
        "test_targets_visible_during_prediction": False,
        "freeze_completed_before_reveal": True,
    }
    (out / "FREEZE_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
