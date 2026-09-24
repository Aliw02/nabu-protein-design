import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd


GENOTYPE_NAMES = {
    "genotype", "genotypes", "mutant", "mutants", "variant", "variants",
    "sequence", "sequences", "aas", "aa", "mutated_sequence", "mutant_sequence",
    "mutation", "mutations"
}
FITNESS_NAMES = {
    "fitness", "score", "scores", "dms_score", "phenotype", "brightness",
    "fluorescence", "activity", "value", "mean", "mean_log"
}
COUNT_NAMES = {
    "mutation_count", "mut_count", "n_mutations", "n_muts", "num_mutations",
    "num_muts", "n_mut", "number_of_mutations"
}
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


def build_mutant_strings(frame, genotype_col, mutation_count):
    notation = frame[genotype_col].map(parse_mutation_notation)
    notation_usable = notation.map(len).gt(0) | mutation_count.fillna(-1).eq(0)

    if notation_usable.mean() > 0.95:
        return notation.map(lambda tokens: ":".join(tokens)), "mutation_notation"

    seq = frame[genotype_col].map(clean_sequence)
    wt_rows = frame[mutation_count.fillna(-1).eq(0)].copy()
    if wt_rows.empty:
        raise RuntimeError(
            "Sequence-mode genotype detected but no mutation_count==0 WT row exists."
        )

    wt = wt_rows[genotype_col].map(clean_sequence)
    wt = wt[wt.str.len() > 0].iloc[0]

    if not seq.map(len).eq(len(wt)).all():
        raise RuntimeError("Inconsistent sequence lengths in eqFP611 dataset.")

    def encode(sequence):
        tokens = []
        for idx, (a, b) in enumerate(zip(wt, sequence), start=1):
            if a != b:
                tokens.append(f"{a}{idx}{b}")
        return ":".join(tokens)

    encoded = seq.map(encode)
    derived = encoded.map(lambda x: 0 if x == "" else len(x.split(":")))
    mismatch = (
        mutation_count.notna()
        & (derived.astype(int) != mutation_count.fillna(-1).astype(int))
    )
    if mismatch.any():
        raise RuntimeError(
            f"Mutation-count mismatch on {int(mismatch.sum())} rows."
        )

    return encoded, "sequence_relative_to_WT"


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
        mutation_count = pd.to_numeric(frame[count_col], errors="coerce").astype("Int64")
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

    visible = frame[frame["bucket"] <= 6].copy()
    hidden = frame[frame["bucket"] >= 7].copy()

    if len(visible) < 500 or len(hidden) < 100:
        raise RuntimeError(
            f"Insufficient 3-5 mutation split: visible={len(visible)}, hidden={len(hidden)}"
        )

    visible_out = pd.DataFrame(
        {
            "candidate_id": visible["candidate_id"],
            "mutant": visible["_mutant"],
            "DMS_score": visible["_fitness"].astype(float),
            "mutation_count": visible["_mutation_count"].astype(int),
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

    visible_path = out / "VISIBLE.csv"
    hidden_ids_path = out / "HIDDEN_IDS.csv"
    hidden_truth_path = out / ".HIDDEN_TRUTH.sealed.csv"

    visible_out.to_csv(visible_path, index=False)
    hidden_ids_out.to_csv(hidden_ids_path, index=False)
    hidden_truth_out.to_csv(hidden_truth_path, index=False)

    manifest = {
        "version": "NABU_V8_3_EQFP611_SEAL_V1",
        "source_file": source.name,
        "source_sha256": sha256_file(source),
        "sheet": sheet,
        "genotype_column": str(genotype_col),
        "fitness_column_name_only": str(fitness_col),
        "mutation_count_source": count_source,
        "genotype_mode": genotype_mode,
        "candidate_scope": "mutation_count in {3,4,5}",
        "split_rule": "FNV1a32(candidate_id) mod 10; visible 0-6, hidden 7-9",
        "visible_rows": int(len(visible_out)),
        "hidden_rows": int(len(hidden_ids_out)),
        "visible_by_mutation_count": {
            str(k): int(v)
            for k, v in visible_out["mutation_count"].value_counts().sort_index().items()
        },
        "hidden_by_mutation_count": {
            str(k): int(v)
            for k, v in hidden_ids_out["mutation_count"].value_counts().sort_index().items()
        },
        "visible_sha256": sha256_file(visible_path),
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
    parser.add_argument("--out", default="eqfp611_sealed_input")
    args = parser.parse_args()
    main(args.xlsx, args.out)
