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


def main(source_csv: str, output_csv: str):
    source = pd.read_csv(source_csv)

    if "variant" not in source.columns or "fitness" not in source.columns:
        raise SystemExit(
            "Expected source columns 'variant' and 'fitness'; "
            f"found {list(source.columns)}"
        )

    data = source[["variant", "fitness"]].dropna().copy()
    data["variant"] = data["variant"].astype(str).str.strip().str.upper()

    valid_length = data["variant"].str.len() == 4
    valid_alphabet = data["variant"].map(
        lambda value: all(aa in STANDARD_AA for aa in value)
    )
    data = data[valid_length & valid_alphabet].copy()

    output = pd.DataFrame(
        {
            "mutant": data["variant"].map(encode_genotype),
            "DMS_score": pd.to_numeric(data["fitness"], errors="coerce"),
            "source_variant": data["variant"],
        }
    ).dropna(subset=["DMS_score"])

    output_path = Path(output_csv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(output_path, index=False)

    print(f"rows={len(output)}")
    print(f"output={output_path}")
    print(output.head(5).to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("source_csv")
    parser.add_argument(
        "--out",
        default="GB1_Wu2016_NABU_V8_1.csv",
    )
    args = parser.parse_args()
    main(args.source_csv, args.out)
