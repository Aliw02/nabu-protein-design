import argparse
from pathlib import Path

import pandas as pd

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def encode_genotype(genotype: str) -> str:
    genotype = str(genotype).strip().upper()
    if len(genotype) != 4:
        raise ValueError(f"Expected four-letter genotype, got: {genotype!r}")
    if any(aa not in STANDARD_AA for aa in genotype):
        raise ValueError(f"Non-standard amino acid in genotype: {genotype!r}")
    return ":".join(
        f"X{position}{aa}"
        for position, aa in enumerate(genotype, start=1)
    )


def main(source_csv, variant_column, fitness_column, output_csv):
    source = pd.read_csv(source_csv)

    if variant_column not in source.columns or fitness_column not in source.columns:
        raise SystemExit(
            f"Expected columns {variant_column!r}/{fitness_column!r}; "
            f"found {list(source.columns)}"
        )

    data = source[[variant_column, fitness_column]].dropna().copy()
    data[variant_column] = (
        data[variant_column].astype(str).str.strip().str.upper()
    )

    valid = (
        data[variant_column].str.len().eq(4)
        & data[variant_column].map(
            lambda value: all(aa in STANDARD_AA for aa in value)
        )
    )
    data = data[valid].copy()

    output = pd.DataFrame(
        {
            "mutant": data[variant_column].map(encode_genotype),
            "DMS_score": pd.to_numeric(
                data[fitness_column],
                errors="coerce",
            ),
            "source_variant": data[variant_column],
        }
    ).dropna(subset=["DMS_score"])

    Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(output_csv, index=False)

    print(f"rows={len(output)}")
    print(f"columns={list(source.columns)}")
    print(f"output={output_csv}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("source_csv")
    parser.add_argument("--variant-column", required=True)
    parser.add_argument("--fitness-column", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    main(
        args.source_csv,
        args.variant_column,
        args.fitness_column,
        args.out,
    )
