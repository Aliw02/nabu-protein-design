from __future__ import annotations

import pandas as pd

from assembly import MutationVocabulary, ReferenceProtein, STANDARD_AA


GB1_REFERENCE_SEQUENCE = (
    "MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE"
)
GB1_MUTABLE_POSITIONS = (39, 40, 41, 54)
GB1_WT_GENOTYPE = "VDGV"


def gb1_reference() -> ReferenceProtein:
    reference = ReferenceProtein(
        name="GB1_3GB1",
        sequence=GB1_REFERENCE_SEQUENCE,
    )
    observed = "".join(
        reference.residue_at(position)
        for position in GB1_MUTABLE_POSITIONS
    )
    if observed != GB1_WT_GENOTYPE:
        raise RuntimeError(
            f"GB1 reference mismatch: expected {GB1_WT_GENOTYPE}, "
            f"found {observed}."
        )
    return reference


def gb1_vocabulary() -> MutationVocabulary:
    reference = gb1_reference()
    tokens = []
    for position in GB1_MUTABLE_POSITIONS:
        source = reference.residue_at(position)
        for target in sorted(STANDARD_AA):
            if target == source:
                continue
            tokens.append(f"{source}{position}{target}")
    return MutationVocabulary.from_tokens(reference, tokens)


def genotype_to_mutation_set(genotype: str) -> tuple[str, ...]:
    genotype = str(genotype).strip().upper()
    if len(genotype) != len(GB1_MUTABLE_POSITIONS):
        raise ValueError(
            f"Expected four-letter GB1 genotype, got {genotype!r}."
        )
    if any(aa not in STANDARD_AA for aa in genotype):
        raise ValueError(f"Non-standard GB1 genotype: {genotype!r}")

    reference = gb1_reference()
    tokens = []
    for position, target in zip(GB1_MUTABLE_POSITIONS, genotype):
        source = reference.residue_at(position)
        if target != source:
            tokens.append(f"{source}{position}{target}")
    return reference.canonicalize_set(tokens)


def adapt_gb1_frame(
    source: pd.DataFrame,
    variant_column: str = "variant",
    fitness_column: str = "fitness",
) -> pd.DataFrame:
    required = {variant_column, fitness_column}
    missing = required - set(source.columns)
    if missing:
        raise ValueError(f"GB1 source missing columns: {sorted(missing)}")

    rows = []
    for row in source[[variant_column, fitness_column]].itertuples(index=False):
        genotype = str(row[0]).strip().upper()
        mutation_set = genotype_to_mutation_set(genotype)
        if len(mutation_set) == 0:
            # The frozen campaign parser has no WT/empty-mutation token.
            # WT is excluded consistently from the 2C candidate universe.
            continue
        label = float(row[1])
        candidate_id = ":".join(mutation_set)
        rows.append(
            {
                "candidate_id": candidate_id,
                "mutant": candidate_id,
                "mutation_set": mutation_set,
                "DMS_score": label,
                "source_variant": genotype,
            }
        )

    frame = pd.DataFrame(rows)
    frame = frame.drop_duplicates("candidate_id", keep="first")
    frame = frame.sort_values("candidate_id").reset_index(drop=True)
    if frame["candidate_id"].duplicated().any():
        raise RuntimeError("GB1 adapter produced duplicate candidate IDs.")
    return frame
