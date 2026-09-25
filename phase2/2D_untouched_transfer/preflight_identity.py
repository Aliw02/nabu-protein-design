from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
from collections import Counter
from itertools import combinations
from pathlib import Path
import zipfile

import pandas as pd

from nabu_protein.higher_order import crossfit_fold

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
AAV_ALLOWED = STANDARD_AA | {"*"}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _summary(values):
    counter = Counter(values)
    return {str(key): int(counter[key]) for key in sorted(counter, key=str)}


def _find_aav_member(names: list[str]) -> str:
    csvs = [name for name in names if name.lower().endswith(".csv")]
    sampled = [name for name in csvs if "sampled" in name.lower()]
    preferred = [
        name for name in sampled
        if "regression" in name.lower()
    ]
    candidates = preferred or sampled
    if len(candidates) != 1:
        raise RuntimeError(
            "Expected exactly one AAV sampled regression CSV; "
            f"candidates={candidates}, all_csvs={csvs}"
        )
    return candidates[0]


def _read_aav_identities(zip_path: Path):
    with zipfile.ZipFile(zip_path) as archive:
        names = archive.namelist()
        member = _find_aav_member(names)
        with archive.open(member) as handle:
            frame = pd.read_csv(
                handle,
                usecols=lambda name: name in {
                    "sequence", "set", "validation"
                },
            )
    required = {"sequence", "set", "validation"}
    if set(frame.columns) != required:
        raise RuntimeError(
            f"AAV identity schema mismatch: {list(frame.columns)}"
        )
    return frame, member, names


def _read_fasta_sequence(path: Path) -> str:
    lines = [
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith(">")
    ]
    sequence = "".join(lines).upper()
    if not sequence:
        raise RuntimeError("Empty AAV reference FASTA.")
    return sequence


def _sequence_profile(frame: pd.DataFrame) -> dict:
    seqs = frame["sequence"].astype(str).str.strip().str.upper()
    alphabet = sorted(set("".join(seqs.tolist())))
    return {
        "rows": int(len(frame)),
        "unique_sequences": int(seqs.nunique()),
        "duplicate_sequence_rows": int(len(frame) - seqs.nunique()),
        "set_counts": {
            str(key): int(value)
            for key, value in frame["set"].astype(str).value_counts().sort_index().items()
        },
        "validation_counts": {
            str(key): int(value)
            for key, value in frame["validation"].astype(str).value_counts().sort_index().items()
        },
        "length_counts": _summary(seqs.str.len().tolist()),
        "alphabet": alphabet,
    }


def _aav_profile(frame: pd.DataFrame, reference: str) -> dict:
    profile = _sequence_profile(frame)

    # Match the released protein-uq loader exactly for the BO candidate pool:
    # set == "train" and validation is NaN. Validation=True rows are held out.
    set_text = frame["set"].astype(str).str.lower()
    published_pool = frame[
        set_text.eq("train") & frame["validation"].isna()
    ].copy()
    seqs = published_pool["sequence"].astype(str).str.strip().str.upper()

    same_length = seqs.str.len().eq(len(reference))
    allowed = seqs.map(lambda seq: set(seq).issubset(AAV_ALLOWED))
    representable = same_length & allowed

    profile.update(
        {
            "reference_length": int(len(reference)),
            "published_bo_pool_rule": "set == train AND validation is NaN",
            "published_bo_pool_rows": int(len(published_pool)),
            "published_bo_pool_length_counts": _summary(seqs.str.len().tolist()),
            "published_bo_pool_alphabet": sorted(set("".join(seqs.tolist()))),
            "representable_rows_frozen_token_grammar": int(representable.sum()),
            "unrepresentable_rows_frozen_token_grammar": int((~representable).sum()),
            "representable_fraction": (
                float(representable.mean()) if len(representable) else 0.0
            ),
            "complete_pool_representable": bool(
                len(representable) > 0 and representable.all()
            ),
            "representation_rule": (
                "published BO pool only; same length as pinned WT and alphabet "
                "limited to 20 standard amino acids plus '*' deletion token"
            ),
        }
    )
    return profile


def _derive_ired_reference(frame: pd.DataFrame) -> dict:
    seqs = frame["sequence"].astype(str).str.strip().str.upper()
    train = frame[frame["set"].astype(str).str.lower().eq("train")].copy()
    train_seqs = train["sequence"].astype(str).str.strip().str.upper()
    lengths = train_seqs.str.len()
    fixed = bool(len(lengths) > 0 and lengths.nunique() == 1)
    standard = bool(train_seqs.map(lambda seq: set(seq).issubset(STANDARD_AA)).all())
    result = {
        "train_fixed_length": fixed,
        "train_standard_amino_acids_only": standard,
        "reference_derivation": "per-position modal residue over official train identities",
        "reference_present_as_exact_row": False,
        "reference_unique": False,
        "derived_reference": None,
    }
    if not fixed or not standard or len(train_seqs) == 0:
        return result

    length = int(lengths.iloc[0])
    reference = "".join(
        Counter(seq[index] for seq in train_seqs).most_common(1)[0][0]
        for index in range(length)
    )
    matches = int((seqs == reference).sum())
    result.update(
        {
            "derived_reference": reference,
            "reference_present_as_exact_row": bool(matches >= 1),
            "reference_exact_row_count": matches,
            "reference_unique": bool(matches == 1),
        }
    )
    return result



def _mutation_set(sequence: str, reference: str) -> tuple[str, ...]:
    if len(sequence) != len(reference):
        raise ValueError("Sequence/reference length mismatch.")
    return tuple(
        f"{source}{index}{target}"
        for index, (source, target) in enumerate(
            zip(reference, sequence),
            start=1,
        )
        if source != target
    )


def _ired_structural_scoreability(
    frame: pd.DataFrame,
    reference: str,
) -> dict:
    set_text = frame["set"].astype(str).str.lower()
    validation_true = frame["validation"].eq(True)

    fit_frame = frame[
        set_text.eq("train") & ~validation_true
    ].copy()
    validation_frame = frame[validation_true].copy()
    test_frame = frame[set_text.eq("test")].copy()

    fit_sets = [
        _mutation_set(str(sequence).strip().upper(), reference)
        for sequence in fit_frame["sequence"]
    ]
    test_sets = [
        _mutation_set(str(sequence).strip().upper(), reference)
        for sequence in test_frame["sequence"]
    ]

    main_support = Counter(
        token
        for mutation_set in fit_sets
        for token in mutation_set
    )
    pair_support = Counter(
        pair
        for mutation_set in fit_sets
        for pair in combinations(mutation_set, 2)
    )

    def scoreable(mutation_set: tuple[str, ...]) -> bool:
        return bool(
            mutation_set
            and all(token in main_support for token in mutation_set)
            and any(pair in pair_support for pair in combinations(mutation_set, 2))
        )

    scoreable_mask = [scoreable(ms) for ms in test_sets]
    missing_main_rows = sum(
        any(token not in main_support for token in ms)
        for ms in test_sets
    )
    no_supported_pair_rows = sum(
        bool(ms)
        and all(
            pair not in pair_support
            for pair in combinations(ms, 2)
        )
        for ms in test_sets
    )

    fold_counts = Counter(
        crossfit_fold(str(sequence).strip().upper())
        for sequence in fit_frame["sequence"]
    )

    return {
        "fit_train_rows": int(len(fit_frame)),
        "validation_rows": int(len(validation_frame)),
        "test_rows": int(len(test_frame)),
        "fit_mutation_order_counts": _summary(
            [len(ms) for ms in fit_sets]
        ),
        "test_mutation_order_counts": _summary(
            [len(ms) for ms in test_sets]
        ),
        "fit_main_identity_entries": int(len(main_support)),
        "fit_pair_identity_entries": int(len(pair_support)),
        "test_structurally_scoreable_rows": int(sum(scoreable_mask)),
        "test_structurally_unscoreable_rows": int(
            len(scoreable_mask) - sum(scoreable_mask)
        ),
        "test_structural_scoreability_fraction": (
            float(sum(scoreable_mask) / len(scoreable_mask))
            if scoreable_mask else 0.0
        ),
        "test_rows_with_unseen_main_identity": int(missing_main_rows),
        "test_rows_without_any_supported_pair_identity": int(
            no_supported_pair_rows
        ),
        "fit_crossfit_fold_counts": {
            str(index): int(fold_counts.get(index, 0))
            for index in range(5)
        },
        "all_five_fit_crossfit_folds_populated": bool(
            set(fold_counts) == set(range(5))
        ),
    }


def preflight(
    aav_zip: Path,
    aav_reference_fasta: Path,
    ired_gz: Path,
) -> dict:
    aav, member, archive_names = _read_aav_identities(aav_zip)
    aav_reference = _read_fasta_sequence(aav_reference_fasta)

    ired = pd.read_csv(
        ired_gz,
        compression="gzip",
        usecols=lambda name: name in {
            "sequence", "set", "validation"
        },
    )
    required = {"sequence", "set", "validation"}
    if set(ired.columns) != required:
        raise RuntimeError(
            f"IRED identity schema mismatch: {list(ired.columns)}"
        )

    ired_profile = _sequence_profile(ired)
    ired_profile.update(_derive_ired_reference(ired))
    if ired_profile["derived_reference"] is not None:
        ired_profile.update(
            _ired_structural_scoreability(
                ired,
                ired_profile["derived_reference"],
            )
        )

    return {
        "version": "NABU_PHASE2D_IDENTITY_PREFLIGHT_V2",
        "target_values_read": False,
        "sources": {
            "aav_splits_zip_sha256": sha256_file(aav_zip),
            "aav_reference_fasta_sha256": sha256_file(aav_reference_fasta),
            "ired_two_to_many_gz_sha256": sha256_file(ired_gz),
        },
        "aav": {
            "selected_archive_member": member,
            "archive_members": archive_names,
            **_aav_profile(aav, aav_reference),
        },
        "ired": ired_profile,
        "gates": {
            "aav_identity_gate_pass": bool(
                _aav_profile(aav, aav_reference)[
                    "complete_pool_representable"
                ]
            ),
            "ired_identity_gate_pass": bool(
                ired_profile["train_fixed_length"]
                and ired_profile["train_standard_amino_acids_only"]
                and ired_profile["reference_present_as_exact_row"]
                and ired_profile.get(
                    "all_five_fit_crossfit_folds_populated",
                    False,
                )
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--aav-zip", required=True)
    parser.add_argument("--aav-reference", required=True)
    parser.add_argument("--ired-gz", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    result = preflight(
        Path(args.aav_zip),
        Path(args.aav_reference),
        Path(args.ired_gz),
    )
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
