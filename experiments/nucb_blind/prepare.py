from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
DEV = HERE.parent / "v9_pair_transfer_dev"
sys.path.insert(0, str(DEV))

from run_v9_pair_transfer_ired import derive_reference, mutation_set


EXPECTED_MD5 = "cdedaadd57f2b148950b357e9f24c238"


def file_hash(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validation_mask(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series.dtype):
        return series.fillna(False).astype(bool)
    text = series.fillna("").astype(str).str.strip().str.lower()
    return text.isin({"true", "1", "yes", "y", "t"})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source_gz")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    source = Path(args.source_gz)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    md5 = file_hash(source, "md5")
    sha256 = file_hash(source, "sha256")
    if md5 != EXPECTED_MD5:
        raise RuntimeError(
            f"NucB MD5 mismatch: expected {EXPECTED_MD5}, got {md5}"
        )

    frame = pd.read_csv(source, compression="gzip")
    required = {"sequence", "target", "set", "validation"}
    missing = required - set(frame.columns)
    if missing:
        raise RuntimeError(
            f"NucB source missing required columns: {sorted(missing)}"
        )

    frame = frame.copy()
    frame["sequence"] = frame["sequence"].astype(str).str.strip().str.upper()
    frame["target"] = pd.to_numeric(frame["target"], errors="raise")
    frame["row_id"] = [f"nucb-{i:09d}" for i in range(len(frame))]

    set_text = frame["set"].astype(str).str.strip().str.lower()
    val = validation_mask(frame["validation"])

    fit = frame[set_text.eq("train") & ~val].copy()
    validation = frame[val].copy()
    test = frame[set_text.eq("test")].copy()

    if len(fit) == 0 or len(test) == 0:
        raise RuntimeError("Empty fit or test partition.")

    reference = derive_reference(fit["sequence"].tolist())
    fit_mutations = [
        mutation_set(sequence, reference)
        for sequence in fit["sequence"]
    ]
    test_mutations = [
        mutation_set(sequence, reference)
        for sequence in test["sequence"]
    ]

    fit_distribution = Counter(len(x) for x in fit_mutations)
    test_distribution = Counter(len(x) for x in test_mutations)

    if any(order > 2 for order in fit_distribution):
        raise RuntimeError(
            f"ADAPTER_SCOPE_FAIL: fit contains mutation order >2: "
            f"{dict(fit_distribution)}"
        )
    if any(order < 3 for order in test_distribution):
        raise RuntimeError(
            f"ADAPTER_SCOPE_FAIL: test contains mutation order <3: "
            f"{dict(test_distribution)}"
        )

    fit[["row_id", "sequence", "target"]].to_csv(
        out / "fit.csv.gz",
        index=False,
        compression="gzip",
    )
    test[["row_id", "sequence"]].to_csv(
        out / "test_masked.csv.gz",
        index=False,
        compression="gzip",
    )
    test[["row_id", "target"]].to_csv(
        out / "test_reveal.csv.gz",
        index=False,
        compression="gzip",
    )

    manifest = {
        "dataset": "FLIP2_NucB_two_to_many",
        "zenodo_record": "18433203",
        "source_md5": md5,
        "source_sha256": sha256,
        "observed_counts": {
            "fit": int(len(fit)),
            "validation": int(len(validation)),
            "test": int(len(test)),
        },
        "reference_length": int(len(reference)),
        "reference_sha256": hashlib.sha256(
            reference.encode("utf-8")
        ).hexdigest(),
        "fit_mutation_order_distribution": {
            str(k): int(v) for k, v in sorted(fit_distribution.items())
        },
        "test_mutation_order_distribution": {
            str(k): int(v) for k, v in sorted(test_distribution.items())
        },
        "validation_targets_used_for_fit": False,
        "test_targets_excluded_from_prediction_artifact": True,
        "scope_contract_pass": True,
    }
    (out / "SOURCE_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
