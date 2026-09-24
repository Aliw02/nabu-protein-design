import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd


RHLAA_SITE_POSITIONS = [74, 101, 143, 148, 173, 176]
RHLAA_WT_STATES = "RAQLSQ"

GENOTYPE_NAMES = [
    "genotype", "genotypes", "mutant", "mutants", "variant", "variants",
    "sequence", "sequences", "aas", "aa", "mutated_sequence", "mutant_sequence",
    "mutation", "mutations",
]
FITNESS_NAMES = [
    "fitness", "combined", "activity", "normalized_activity",
    "score", "scores", "dms_score", "phenotype", "value", "mean", "mean_log",
]
COUNT_NAMES = [
    "n_mut", "mutation_count", "mut_count", "n_mutations", "n_muts",
    "num_mutations", "num_muts", "number_of_mutations",
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


def find_data_sheet(xlsx_path):
    book = pd.ExcelFile(xlsx_path)
    observed = []
    for sheet in book.sheet_names:
        frame = pd.read_excel(xlsx_path, sheet_name=sheet)
        genotype_col = choose_column(frame.columns, GENOTYPE_NAMES)
        fitness_col = choose_column(frame.columns, FITNESS_NAMES)
        count_col = choose_column(frame.columns, COUNT_NAMES)
        if genotype_col is not None and fitness_col is not None:
            return sheet, frame, genotype_col, fitness_col, count_col
        observed.append(
            {"sheet": sheet, "columns": [str(c) for c in frame.columns]}
        )
    raise RuntimeError(f"Could not identify required columns. Observed: {observed}")


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


def encode_compact_rhla(sequence):
    if len(sequence) != len(RHLAA_SITE_POSITIONS):
        raise ValueError(sequence)
    tokens = []
    for wt, aa, position in zip(
        RHLAA_WT_STATES,
        sequence,
        RHLAA_SITE_POSITIONS,
    ):
        if aa != wt:
            tokens.append(f"{wt}{position}{aa}")
    return ":".join(tokens)


def build_mutant_strings(frame, genotype_col, mutation_count):
    notation = frame[genotype_col].map(parse_mutation_notation)
    notation_usable = notation.map(len).gt(0) | mutation_count.fillna(-1).eq(0)

    if notation_usable.mean() > 0.95:
        encoded = notation.map(lambda tokens: ":".join(tokens))
        return encoded, "mutation_notation"

    sequence = frame[genotype_col].map(clean_sequence)
    lengths = sequence.map(len)

    if lengths.eq(len(RHLAA_SITE_POSITIONS)).mean() > 0.95:
        encoded = sequence.map(encode_compact_rhla)
        derived = encoded.map(lambda x: 0 if x == "" else len(x.split(":")))
        mismatch = (
            mutation_count.notna()
            & (derived.astype(int) != mutation_count.fillna(-1).astype(int))
        )
        if mismatch.any():
            raise RuntimeError(
                f"Compact RhlA genotype mutation-count mismatch on "
                f"{int(mismatch.sum())} rows."
            )
        return encoded, "compact_6_site_rhla"

    wt_rows = frame[mutation_count.fillna(-1).eq(0)].copy()
    if wt_rows.empty:
        raise RuntimeError(
            "Full-sequence genotype detected but no mutation_count==0 WT row exists."
        )

    wt = wt_rows[genotype_col].map(clean_sequence)
    wt = wt[wt.str.len() > 0].iloc[0]
    if not lengths.eq(len(wt)).all():
        raise RuntimeError("Inconsistent full-sequence lengths in RhlA dataset.")

    def encode_full(seq):
        tokens = []
        for idx, (a, b) in enumerate(zip(wt, seq), start=1):
            if a != b:
                tokens.append(f"{a}{idx}{b}")
        return ":".join(tokens)

    encoded = sequence.map(encode_full)
    derived = encoded.map(lambda x: 0 if x == "" else len(x.split(":")))
    mismatch = (
        mutation_count.notna()
        & (derived.astype(int) != mutation_count.fillna(-1).astype(int))
    )
    if mismatch.any():
        raise RuntimeError(
            f"Full-sequence mutation-count mismatch on {int(mismatch.sum())} rows."
        )

    return encoded, "full_sequence_relative_to_WT"


def main(xlsx_path, out_dir):
    source = Path(xlsx_path)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    sheet, frame, genotype_col, fitness_col, count_col = find_data_sheet(source)

    if count_col is None:
        notation = frame[genotype_col].map(parse_mutation_notation)
        if not notation.map(len).gt(0).any():
            raise RuntimeError("No mutation-count column and no mutation notation.")
        mutation_count = notation.map(len).astype("Int64")
        count_source = "derived_from_mutation_notation"
    else:
        mutation_count = pd.to_numeric(
            frame[count_col], errors="coerce"
        ).astype("Int64")
        count_source = str(count_col)

    frame = frame.copy()
    frame["_mutation_count"] = mutation_count
    frame = frame[frame["_mutation_count"].notna()].copy()

    mutant, genotype_mode = build_mutant_strings(
        frame,
        genotype_col,
        frame["_mutation_count"],
    )
    frame["_mutant"] = mutant
    frame["_fitness"] = pd.to_numeric(frame[fitness_col], errors="coerce")
    frame = frame[frame["_fitness"].notna()].copy()

    frame = frame[frame["_mutation_count"].between(3, 5)].copy()
    frame["candidate_id"] = frame["_mutant"].astype(str)
    frame = frame.drop_duplicates("candidate_id", keep="first").copy()
    frame["bucket"] = frame["candidate_id"].map(lambda x: fnv1a32(x) % 10)

    training_pool = frame[frame["bucket"] <= 6].copy()
    hidden = frame[frame["bucket"] >= 7].copy()

    if len(training_pool) < 200 or len(hidden) < 50:
        raise RuntimeError(
            "Insufficient 3-5 mutation split: "
            f"training_pool={len(training_pool)}, hidden={len(hidden)}"
        )

    training_pool_out = pd.DataFrame(
        {
            "candidate_id": training_pool["candidate_id"],
            "mutant": training_pool["_mutant"],
            "DMS_score": training_pool["_fitness"].astype(float),
            "mutation_count": training_pool["_mutation_count"].astype(int),
        }
    )
    hidden_ids_out = pd.DataFrame(
        {
            "candidate_id": hidden["candidate_id"],
            "mutant": hidden["_mutant"],
            "mutation_count": hidden["_mutation_count"].astype(int),
        }
    )
    hidden_truth_out = pd.DataFrame(
        {
            "candidate_id": hidden["candidate_id"],
            "DMS_score": hidden["_fitness"].astype(float),
            "mutation_count": hidden["_mutation_count"].astype(int),
        }
    )

    training_pool_path = out / "TRAINING_POOL.csv"
    hidden_ids_path = out / "HIDDEN_IDS.csv"
    hidden_truth_path = out / ".HIDDEN_TRUTH.sealed.csv"

    training_pool_out.to_csv(training_pool_path, index=False)
    hidden_ids_out.to_csv(hidden_ids_path, index=False)
    hidden_truth_out.to_csv(hidden_truth_path, index=False)

    manifest = {
        "version": "NABU_V8_3_RHLA_SEAL_V1",
        "source_file": source.name,
        "source_sha256": sha256_file(source),
        "sheet": sheet,
        "genotype_column": str(genotype_col),
        "fitness_column_name_only": str(fitness_col),
        "mutation_count_source": count_source,
        "genotype_mode": genotype_mode,
        "candidate_scope": "mutation_count in {3,4,5}",
        "split_rule": "FNV1a32(candidate_id) mod 10; training_pool 0-6, hidden 7-9",
        "training_pool_rows": int(len(training_pool_out)),
        "hidden_rows": int(len(hidden_ids_out)),
        "total_eligible_rows": int(len(training_pool_out) + len(hidden_ids_out)),
        "training_pool_by_mutation_count": {
            str(k): int(v)
            for k, v in training_pool_out[
                "mutation_count"
            ].value_counts().sort_index().items()
        },
        "hidden_by_mutation_count": {
            str(k): int(v)
            for k, v in hidden_ids_out[
                "mutation_count"
            ].value_counts().sort_index().items()
        },
        "training_pool_sha256": sha256_file(training_pool_path),
        "hidden_ids_sha256": sha256_file(hidden_ids_path),
        "hidden_truth_sha256": sha256_file(hidden_truth_path),
        "hidden_truth_values_printed": False,
    }
    (out / "SEAL_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )

    print(json.dumps(manifest, indent=2))
    print("SEALED: hidden fitness values were not printed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("xlsx")
    parser.add_argument("--out", default="rhla_sealed_input")
    args = parser.parse_args()
    main(args.xlsx, args.out)
