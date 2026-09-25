from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from itertools import combinations
from pathlib import Path

import pandas as pd

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def derive_reference(train_sequences: list[str]) -> str:
    if not train_sequences:
        raise ValueError("Empty IRED fit/train identity set.")
    lengths = {len(seq) for seq in train_sequences}
    if len(lengths) != 1:
        raise ValueError(f"IRED fit/train sequences are not fixed-length: {lengths}")
    if any(not set(seq).issubset(STANDARD_AA) for seq in train_sequences):
        raise ValueError("IRED fit/train contains non-standard amino acids.")
    length = next(iter(lengths))
    return "".join(
        Counter(seq[i] for seq in train_sequences).most_common(1)[0][0]
        for i in range(length)
    )


def mutation_set(sequence: str, reference: str) -> tuple[str, ...]:
    sequence = str(sequence).strip().upper()
    if len(sequence) != len(reference):
        raise ValueError("Sequence/reference length mismatch.")
    if not set(sequence).issubset(STANDARD_AA):
        raise ValueError("Non-standard amino acid in IRED sequence.")
    return tuple(
        f"{source}{i + 1}{target}"
        for i, (source, target) in enumerate(zip(reference, sequence))
        if source != target
    )


def scoreability_preflight(source_gz: Path) -> dict:
    frame = pd.read_csv(
        source_gz,
        compression="gzip",
        usecols=lambda name: name in {"sequence", "set", "validation"},
    )
    required = {"sequence", "set", "validation"}
    if set(frame.columns) != required:
        raise RuntimeError(f"IRED identity schema mismatch: {list(frame.columns)}")

    seq = frame["sequence"].astype(str).str.strip().str.upper()
    set_text = frame["set"].astype(str).str.lower()
    validation_true = frame["validation"].fillna(False).astype(bool)

    fit_mask = set_text.eq("train") & ~validation_true
    val_mask = validation_true
    test_mask = set_text.eq("test")

    fit_sequences = seq[fit_mask].tolist()
    reference = derive_reference(fit_sequences)

    fit_sets = [mutation_set(value, reference) for value in fit_sequences]
    test_sequences = seq[test_mask].tolist()
    test_sets = [mutation_set(value, reference) for value in test_sequences]

    main_support = Counter(
        token
        for mutations in fit_sets
        for token in mutations
    )
    pair_support = Counter(
        pair
        for mutations in fit_sets
        for pair in combinations(mutations, 2)
    )

    rows = []
    for sequence, mutations in zip(test_sequences, test_sets):
        missing_main = tuple(token for token in mutations if token not in main_support)
        supported_pairs = sum(
            pair in pair_support
            for pair in combinations(mutations, 2)
        )
        all_main = len(missing_main) == 0
        any_pair = supported_pairs > 0
        scoreable = all_main and any_pair
        rows.append(
            {
                "sequence_sha256": sha256_text(sequence),
                "mutation_count": len(mutations),
                "all_main_supported": all_main,
                "supported_pair_count": int(supported_pairs),
                "scoreable": scoreable,
                "missing_main_count": len(missing_main),
            }
        )

    out = pd.DataFrame(rows)
    mutation_counts = Counter(out["mutation_count"].astype(int).tolist())
    unscoreable = out[~out["scoreable"]]

    return {
        "version": "NABU_PHASE2D_IRED_SCOREABILITY_PREFLIGHT_V1",
        "target_values_read": False,
        "fit_rule": "set == train AND validation != True",
        "validation_rule": "validation == True held out",
        "test_rule": "set == test",
        "fit_count": int(fit_mask.sum()),
        "validation_count": int(val_mask.sum()),
        "test_count": int(test_mask.sum()),
        "reference_sha256": sha256_text(reference),
        "reference_length": len(reference),
        "fit_main_identity_entries": len(main_support),
        "fit_pair_identity_entries": len(pair_support),
        "test_mutation_count_distribution": {
            str(k): int(v) for k, v in sorted(mutation_counts.items())
        },
        "all_main_supported_count": int(out["all_main_supported"].sum()),
        "any_pair_supported_count": int((out["supported_pair_count"] > 0).sum()),
        "scoreable_count": int(out["scoreable"].sum()),
        "unscoreable_count": int((~out["scoreable"]).sum()),
        "scoreable_fraction": float(out["scoreable"].mean()),
        "full_test_scoreable": bool(out["scoreable"].all()),
        "unscoreable_examples": unscoreable.head(20).to_dict(orient="records"),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ired-gz", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    result = scoreability_preflight(Path(args.ired_gz))
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
