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
    "sequence", "sequences", "aas", "aa", "mutated_sequence",
    "mutant_sequence", "mutation", "mutations",
]
COUNT_NAMES = [
    "num_mut", "n_mut", "mutation_count", "mut_count",
    "n_mutations", "n_muts", "num_mutations", "num_muts",
    "number_of_mutations",
]
FITNESS_EXACT = [
    "fitness", "dms_score", "score", "binding", "affinity",
    "h1_binding", "h1_affinity", "neg_log_kd", "minus_log_kd",
    "log_kd", "logkd", "kd",
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


def choose_fitness_column(frame, genotype_col, count_col):
    normalized = {norm(c): c for c in frame.columns}

    for name in FITNESS_EXACT:
        if name in normalized:
            return normalized[name]

    preferred = []
    for key, original in normalized.items():
        if original in {genotype_col, count_col}:
            continue
        if any(token in key for token in ["h1", "affinity", "binding", "fitness"]):
            numeric = pd.to_numeric(frame[original], errors="coerce")
            if numeric.notna().mean() > 0.90:
                preferred.append(original)

    preferred = list(dict.fromkeys(preferred))
    if len(preferred) == 1:
        return preferred[0]

    numeric_candidates = []
    for original in frame.columns:
        if original in {genotype_col, count_col}:
            continue
        key = norm(original)
        if key in {"index", "id", "row", "unnamed_0"}:
            continue
        numeric = pd.to_numeric(frame[original], errors="coerce")
        if numeric.notna().mean() > 0.95:
            numeric_candidates.append(original)

    if len(numeric_candidates) == 1:
        return numeric_candidates[0]

    raise RuntimeError(
        "Could not identify one unambiguous CR9114-H1 fitness column. "
        f"preferred={preferred}; numeric_candidates={numeric_candidates}; "
        f"columns={[str(c) for c in frame.columns]}"
    )


def find_data_sheet(xlsx_path):
    book = pd.ExcelFile(xlsx_path)
    observed = []

    for sheet in book.sheet_names:
        frame = pd.read_excel(xlsx_path, sheet_name=sheet)
        genotype_col = choose_column(frame.columns, GENOTYPE_NAMES)
        count_col = choose_column(frame.columns, COUNT_NAMES)

        fitness_col = None
        if genotype_col is not None and count_col is not None:
            try:
                fitness_col = choose_fitness_column(
                    frame,
                    genotype_col,
                    count_col,
                )
            except RuntimeError:
                fitness_col = None

        if (
            genotype_col is not None
            and count_col is not None
            and fitness_col is not None
        ):
            return sheet, frame, genotype_col, fitness_col, count_col

        observed.append({
            "sheet": sheet,
            "columns": [str(c) for c in frame.columns],
            "genotype_column": None if genotype_col is None else str(genotype_col),
            "count_column": None if count_col is None else str(count_col),
            "fitness_column": None if fitness_col is None else str(fitness_col),
        })

    raise RuntimeError(
        "Could not find one CR9114-H1 sheet with genotype, mutation count, "
        f"and unambiguous H1 binding score. Observed={observed}"
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
                "Mutation-notation counts disagree with num_mut on too many rows: "
                f"{int(mismatch.sum())}/{len(frame)}"
            )
        return encoded, "mutation_notation", None

    sequences = frame[genotype_col].map(clean_sequence)
    if sequences.map(len).eq(0).any():
        raise RuntimeError("Some CR9114 genotype identities are empty after cleaning.")

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
        raise RuntimeError("Reference length does not match CR9114 genotype identities.")

    def encode(sequence):
        tokens = []
        for idx, (ref_state, state) in enumerate(zip(reference, sequence), start=1):
            if ref_state != state:
                tokens.append(f"{ref_state}{idx}{state}")
        return ":".join(tokens)

    encoded = sequences.map(encode)
    derived = encoded.map(
        lambda x: 0 if x == "" else len(x.split(":"))
    ).astype(int)

    mismatch = (
        mutation_count.notna()
        & (derived != mutation_count.fillna(-1).astype(int))
    )

    if float(mismatch.mean()) > 0.01:
        raise RuntimeError(
            "Identity-derived mutation counts disagree with num_mut on too many rows. "
            f"mismatches={int(mismatch.sum())}/{len(frame)}; "
            f"reference_width={len(reference)}; "
            f"reference_source={reference_source}"
        )

    mode = (
        "sequence_or_compact_identity_relative_to_reference:"
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
        frame["_mutation_count"].notna()
        & frame["_fitness"].notna()
    ].copy()

    mutant, genotype_mode, identity_reference = build_mutant_strings(
        frame,
        genotype_col,
        frame["_mutation_count"],
    )
    frame["_mutant"] = mutant

    frame = frame[frame["_mutation_count"].between(1, 5)].copy()
    frame["candidate_id"] = frame["_mutant"].astype(str)
    frame = frame[
        frame["candidate_id"].str.len().gt(0)
    ].drop_duplicates("candidate_id", keep="first").copy()

    frame["bucket"] = frame["candidate_id"].map(
        lambda x: fnv1a32(x) % 10
    )

    training_pool = frame[frame["bucket"] <= 6].copy()

    hidden_test = frame[
        (frame["bucket"] >= 7)
        & frame["_mutation_count"].isin([4, 5])
    ].copy()

    unused_holdout = frame[
        (frame["bucket"] >= 7)
        & frame["_mutation_count"].isin([1, 2, 3])
    ].copy()

    if len(training_pool) < 1000:
        raise RuntimeError(
            f"CR9114 training pool unexpectedly small: {len(training_pool)}"
        )
    if len(hidden_test) < 500:
        raise RuntimeError(
            f"CR9114 hidden 4/5-mutant test unexpectedly small: {len(hidden_test)}"
        )

    training_out = pd.DataFrame({
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

    training_path = out / "TRAINING_POOL.csv"
    hidden_ids_path = out / "HIDDEN_4_5_IDS.csv"
    hidden_truth_path = out / ".HIDDEN_4_5_TRUTH.sealed.csv"

    training_out.to_csv(training_path, index=False)
    hidden_ids_out.to_csv(hidden_ids_path, index=False)
    hidden_truth_out.to_csv(hidden_truth_path, index=False)

    total_eligible_rows = int(len(frame))

    manifest = {
        "version": "NABU_V8_3_CR9114_H1_FINAL_GATE_SEAL_V1",
        "source_file": source.name,
        "source_sha256": sha256_file(source),
        "sheet": sheet,
        "genotype_column": str(genotype_col),
        "fitness_column_name_only": str(fitness_col),
        "mutation_count_column": str(count_col),
        "genotype_mode": genotype_mode,
        "identity_reference_if_inferred": identity_reference,
        "candidate_scope": "mutation_count in {1,2,3,4,5}",
        "hidden_evaluation_scope": "mutation_count in {4,5}",
        "split_rule": (
            "FNV1a32(candidate_id) mod 10; training 0-6; "
            "hidden 4/5-mutants 7-9"
        ),
        "total_eligible_rows": total_eligible_rows,
        "training_pool_rows": int(len(training_out)),
        "hidden_4_5_rows": int(len(hidden_ids_out)),
        "unused_lower_order_holdout_rows": int(len(unused_holdout)),
        "training_pool_by_mutation_count": {
            str(k): int(v)
            for k, v in training_out[
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
        "training_pool_sha256": sha256_file(training_path),
        "hidden_ids_sha256": sha256_file(hidden_ids_path),
        "hidden_truth_sha256": sha256_file(hidden_truth_path),
        "hidden_truth_values_printed": False,
    }

    if max(manifest["budget_target_counts"].values()) > len(training_out):
        raise RuntimeError(
            "40% budget exceeds deterministic CR9114 training-pool capacity."
        )

    (out / "SEAL_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )

    print(json.dumps(manifest, indent=2))
    print("SEALED: hidden CR9114 4/5-mutant fitness values were not printed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("xlsx")
    parser.add_argument(
        "--out",
        default="cr9114_h1_final_gate_sealed_input",
    )
    args = parser.parse_args()
    main(args.xlsx, args.out)
