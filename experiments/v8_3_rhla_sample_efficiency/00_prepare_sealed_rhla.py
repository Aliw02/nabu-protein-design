import argparse
import hashlib
import json
import math
import re
from collections import Counter
from pathlib import Path

import pandas as pd


GENOTYPE_NAMES = [
    "genotype", "genotypes", "mutant", "mutants", "variant", "variants",
    "sequence", "sequences", "aas", "aa", "mutated_sequence", "mutant_sequence",
    "mutation", "mutations",
]
COUNT_NAMES = [
    "n_mut", "mutation_count", "mut_count", "n_mutations", "n_muts",
    "num_mutations", "num_muts", "number_of_mutations",
]
ACTIVITY_EXACT = [
    "activity",
    "normalized_activity",
    "overall_activity",
    "total_activity",
    "enzyme_activity",
]
TOKEN_RE = re.compile(r"[A-Z][0-9]+[A-Z*]")


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def fnv1a32(text):
    h = 2166136261
    for byte in str(text).encode("utf-8"):
        h ^= byte
        h = (h * 16777619) & 0xFFFFFFFF
    return h


def norm(name):
    return re.sub(r"[^a-z0-9]+", "_", str(name).strip().lower()).strip("_")


def choose_column(columns, accepted):
    normalized = {norm(c): c for c in columns}
    for name in accepted:
        if name in normalized:
            return normalized[name]
    return None


def choose_activity_column(columns):
    normalized = {norm(c): c for c in columns}

    for name in ACTIVITY_EXACT:
        if name in normalized:
            return normalized[name]

    candidates = []
    for normalized_name, original in normalized.items():
        if "activity" not in normalized_name:
            continue
        if "specific" in normalized_name:
            continue
        if "select" in normalized_name:
            continue
        candidates.append(original)

    if len(candidates) == 1:
        return candidates[0]

    raise RuntimeError(
        "Could not identify one unambiguous RhlA overall-activity column. "
        f"Activity-like candidates={candidates}; all columns={list(columns)}"
    )


def find_data_sheet(xlsx_path):
    book = pd.ExcelFile(xlsx_path)
    observed = []

    for sheet in book.sheet_names:
        frame = pd.read_excel(xlsx_path, sheet_name=sheet)
        genotype_col = choose_column(frame.columns, GENOTYPE_NAMES)
        count_col = choose_column(frame.columns, COUNT_NAMES)

        try:
            fitness_col = choose_activity_column(frame.columns)
        except RuntimeError:
            fitness_col = None

        if genotype_col is not None and count_col is not None and fitness_col is not None:
            return sheet, frame, genotype_col, fitness_col, count_col

        observed.append({
            "sheet": sheet,
            "columns": [str(c) for c in frame.columns],
            "genotype_column": None if genotype_col is None else str(genotype_col),
            "count_column": None if count_col is None else str(count_col),
            "activity_column": None if fitness_col is None else str(fitness_col),
        })

    raise RuntimeError(
        "Could not find a sheet with genotype, mutation count, and one "
        f"unambiguous overall-activity column. Observed={observed}"
    )


def parse_mutation_notation(value):
    tokens = TOKEN_RE.findall(str(value).upper())
    return tuple(
        sorted(
            set(tokens),
            key=lambda token: (int(token[1:-1]), token),
        )
    )


def clean_sequence(value):
    return re.sub(r"[^A-Za-z*]", "", str(value)).upper()


def infer_reference_from_identities(sequences):
    lengths = sequences.map(len)
    if lengths.nunique() != 1:
        raise RuntimeError(
            "Sequence-style genotypes have inconsistent lengths: "
            f"{sorted(lengths.unique().tolist())}"
        )

    width = int(lengths.iloc[0])
    if width < 1:
        raise RuntimeError("Empty sequence-style genotype.")

    chars = []
    for idx in range(width):
        counts = Counter(seq[idx] for seq in sequences)
        state, _ = counts.most_common(1)[0]
        chars.append(state)
    return "".join(chars)


def build_mutant_strings(frame, genotype_col, mutation_count):
    notation = frame[genotype_col].map(parse_mutation_notation)
    notation_usable = notation.map(len).gt(0)

    if notation_usable.mean() > 0.95:
        encoded = notation.map(lambda tokens: ":".join(tokens))
        derived = notation.map(len).astype(int)
        mismatch = (
            mutation_count.notna()
            & (derived != mutation_count.fillna(-1).astype(int))
        )
        if mismatch.mean() > 0.01:
            raise RuntimeError(
                "Mutation-notation counts disagree with n_mut on too many rows: "
                f"{int(mismatch.sum())}/{len(frame)}"
            )
        return encoded, "mutation_notation", None

    sequences = frame[genotype_col].map(clean_sequence)
    if sequences.map(len).eq(0).any():
        raise RuntimeError("Some sequence-style genotypes are empty after cleaning.")

    zero_rows = frame[mutation_count.fillna(-1).eq(0)].copy()
    reference = None
    reference_source = None

    if not zero_rows.empty:
        zero_sequences = zero_rows[genotype_col].map(clean_sequence)
        zero_sequences = zero_sequences[zero_sequences.map(len) > 0]
        if not zero_sequences.empty:
            reference = zero_sequences.iloc[0]
            reference_source = "explicit_zero_mutant"

    if reference is None:
        reference = infer_reference_from_identities(sequences)
        reference_source = "identity_only_position_mode"

    if not sequences.map(len).eq(len(reference)).all():
        raise RuntimeError("Reference length does not match all genotype identities.")

    def encode(sequence):
        tokens = []
        for idx, (ref_state, state) in enumerate(zip(reference, sequence), start=1):
            if ref_state != state:
                tokens.append(f"{ref_state}{idx}{state}")
        return ":".join(tokens)

    encoded = sequences.map(encode)
    derived = encoded.map(lambda x: 0 if x == "" else len(x.split(":"))).astype(int)
    mismatch = (
        mutation_count.notna()
        & (derived != mutation_count.fillna(-1).astype(int))
    )

    mismatch_fraction = float(mismatch.mean())
    if mismatch_fraction > 0.01:
        raise RuntimeError(
            "Identity-derived mutation counts disagree with n_mut on too many rows. "
            f"mismatches={int(mismatch.sum())}/{len(frame)} "
            f"({mismatch_fraction:.3%}); genotype width={len(reference)}; "
            f"reference_source={reference_source}"
        )

    mode = (
        "full_or_compact_sequence_relative_to_identity_reference:"
        + reference_source
    )
    return encoded, mode, reference


def main(xlsx_path, out_dir):
    source = Path(xlsx_path)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    sheet, frame, genotype_col, fitness_col, count_col = find_data_sheet(source)

    frame = frame.copy()
    frame["_mutation_count"] = pd.to_numeric(
        frame[count_col],
        errors="coerce",
    ).astype("Int64")
    frame["_fitness"] = pd.to_numeric(
        frame[fitness_col],
        errors="coerce",
    )
    frame = frame[
        frame["_mutation_count"].notna() & frame["_fitness"].notna()
    ].copy()

    mutant, genotype_mode, identity_reference = build_mutant_strings(
        frame,
        genotype_col,
        frame["_mutation_count"],
    )
    frame["_mutant"] = mutant

    frame = frame[frame["_mutation_count"].between(1, 3)].copy()
    frame["candidate_id"] = frame["_mutant"].astype(str)
    frame = frame[
        frame["candidate_id"].str.len().gt(0)
    ].copy()
    frame = frame.drop_duplicates("candidate_id", keep="first").copy()

    frame["bucket"] = frame["candidate_id"].map(
        lambda x: fnv1a32(x) % 10
    )

    training_pool = frame[frame["bucket"] <= 6].copy()
    hidden_test = frame[
        (frame["bucket"] >= 7)
        & frame["_mutation_count"].eq(3)
    ].copy()
    unused_lower_order_holdout = frame[
        (frame["bucket"] >= 7)
        & frame["_mutation_count"].lt(3)
    ].copy()

    if len(training_pool) < 100:
        raise RuntimeError(
            f"Training pool too small: {len(training_pool)}"
        )
    if len(hidden_test) < 50:
        raise RuntimeError(
            f"Hidden triple-mutant test too small: {len(hidden_test)}"
        )

    training_pool_out = pd.DataFrame({
        "candidate_id": training_pool["candidate_id"],
        "mutant": training_pool["_mutant"],
        "DMS_score": training_pool["_fitness"].astype(float),
        "mutation_count": training_pool["_mutation_count"].astype(int),
    })

    hidden_ids_out = pd.DataFrame({
        "candidate_id": hidden_test["candidate_id"],
        "mutant": hidden_test["_mutant"],
        "mutation_count": hidden_test["_mutation_count"].astype(int),
    })

    hidden_truth_out = pd.DataFrame({
        "candidate_id": hidden_test["candidate_id"],
        "DMS_score": hidden_test["_fitness"].astype(float),
        "mutation_count": hidden_test["_mutation_count"].astype(int),
    })

    training_pool_path = out / "TRAINING_POOL.csv"
    hidden_ids_path = out / "HIDDEN_TRIPLE_IDS.csv"
    hidden_truth_path = out / ".HIDDEN_TRIPLE_TRUTH.sealed.csv"

    training_pool_out.to_csv(training_pool_path, index=False)
    hidden_ids_out.to_csv(hidden_ids_path, index=False)
    hidden_truth_out.to_csv(hidden_truth_path, index=False)

    total_eligible_rows = int(len(frame))

    manifest = {
        "version": "NABU_V8_3_RHLA_SAMPLE_EFFICIENCY_SEAL_V1",
        "source_file": source.name,
        "source_sha256": sha256_file(source),
        "sheet": sheet,
        "genotype_column": str(genotype_col),
        "fitness_objective": "overall enzyme activity",
        "fitness_column_name_only": str(fitness_col),
        "mutation_count_column": str(count_col),
        "genotype_mode": genotype_mode,
        "identity_reference_if_inferred": identity_reference,
        "candidate_scope": "mutation_count in {1,2,3}",
        "hidden_evaluation_scope": "triple mutants only",
        "split_rule": (
            "FNV1a32(candidate_id) mod 10; training pool 0-6; "
            "hidden triple test 7-9"
        ),
        "total_eligible_rows": total_eligible_rows,
        "training_pool_rows": int(len(training_pool_out)),
        "hidden_triple_rows": int(len(hidden_ids_out)),
        "unused_lower_order_holdout_rows": int(len(unused_lower_order_holdout)),
        "training_pool_by_mutation_count": {
            str(k): int(v)
            for k, v in training_pool_out[
                "mutation_count"
            ].value_counts().sort_index().items()
        },
        "hidden_test_by_mutation_count": {
            str(k): int(v)
            for k, v in hidden_ids_out[
                "mutation_count"
            ].value_counts().sort_index().items()
        },
        "budget_target_counts": {
            str(p): int(math.ceil(total_eligible_rows * p / 100.0))
            for p in [5, 10, 20, 40]
        },
        "training_pool_sha256": sha256_file(training_pool_path),
        "hidden_ids_sha256": sha256_file(hidden_ids_path),
        "hidden_truth_sha256": sha256_file(hidden_truth_path),
        "hidden_truth_values_printed": False,
    }

    max_target = max(manifest["budget_target_counts"].values())
    if max_target > len(training_pool_out):
        raise RuntimeError(
            "40% budget exceeds deterministic training-pool capacity: "
            f"target={max_target}, pool={len(training_pool_out)}"
        )

    (out / "SEAL_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )

    print(json.dumps(manifest, indent=2))
    print("SEALED: hidden triple fitness values were not printed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("xlsx")
    parser.add_argument(
        "--out",
        default="rhla_sample_efficiency_sealed_input",
    )
    args = parser.parse_args()
    main(args.xlsx, args.out)
