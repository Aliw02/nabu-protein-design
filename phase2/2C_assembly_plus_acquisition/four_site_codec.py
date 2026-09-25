from __future__ import annotations

from dataclasses import dataclass
import hashlib

import pandas as pd

from nabu_protein.higher_order import crossfit_fold, fnv1a32
from nabu_protein.v83 import canonicalize_mutation_set


STANDARD_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
POOL_NAMESPACE = "NABU_PHASE2C_ACQUISITION_POOL_V1"
BOOTSTRAP_NAMESPACE = "NABU_PHASE2C_BOOTSTRAP_V1"


@dataclass(frozen=True)
class FourSiteCodec:
    reference_genotype: str

    def __post_init__(self):
        reference = str(self.reference_genotype).strip().upper()
        if len(reference) != 4:
            raise ValueError(
                "reference_genotype must contain exactly four residues."
            )
        invalid = sorted(set(reference) - STANDARD_AA)
        if invalid:
            raise ValueError(
                "reference_genotype contains non-standard residues: "
                f"{invalid}"
            )
        object.__setattr__(self, "reference_genotype", reference)

    def encode(self, genotype: str) -> tuple[str, ...]:
        genotype = str(genotype).strip().upper()
        if len(genotype) != 4:
            raise ValueError(
                f"Expected four-residue genotype, got {genotype!r}."
            )
        invalid = sorted(set(genotype) - STANDARD_AA)
        if invalid:
            raise ValueError(
                f"Genotype contains non-standard residues: {invalid}"
            )

        tokens = []
        for position, (source, target) in enumerate(
            zip(self.reference_genotype, genotype),
            start=1,
        ):
            if source == target:
                continue
            tokens.append(f"{source}{position}{target}")
        return canonicalize_mutation_set(tokens)

    def decode(self, mutation_set) -> str:
        residues = list(self.reference_genotype)
        canonical = canonicalize_mutation_set(mutation_set)
        used_positions = set()

        for token in canonical:
            source = token[0]
            position = int(token[1:-1])
            target = token[-1]

            if position < 1 or position > 4:
                raise ValueError(
                    f"Mutation {token} is outside the four-site benchmark."
                )
            if position in used_positions:
                raise ValueError(
                    f"Multiple substitutions target position {position}."
                )
            used_positions.add(position)

            reference_source = self.reference_genotype[position - 1]
            if source != reference_source:
                raise ValueError(
                    f"Mutation {token} source {source} does not match "
                    f"reference {reference_source} at position {position}."
                )
            residues[position - 1] = target

        return "".join(residues)


def load_four_site_landscape(
    path: str,
    variant_column: str,
    fitness_column: str,
    reference_genotype: str,
) -> tuple[pd.DataFrame, FourSiteCodec]:
    source = pd.read_csv(path)
    required = {variant_column, fitness_column}
    missing = required - set(source.columns)
    if missing:
        raise ValueError(
            f"Landscape {path} missing columns: {sorted(missing)}"
        )

    codec = FourSiteCodec(reference_genotype)

    frame = source[[variant_column, fitness_column]].copy()
    frame.columns = ["candidate_id", "DMS_score"]
    frame["candidate_id"] = (
        frame["candidate_id"].astype(str).str.strip().str.upper()
    )
    frame["DMS_score"] = pd.to_numeric(
        frame["DMS_score"],
        errors="coerce",
    )

    valid_identity = (
        frame["candidate_id"].str.len().eq(4)
        & frame["candidate_id"].map(
            lambda value: all(
                residue in STANDARD_AA
                for residue in value
            )
        )
    )
    frame = frame[
        valid_identity
    ].dropna(subset=["DMS_score"]).copy()

    if frame["candidate_id"].duplicated().any():
        frame = (
            frame.groupby(
                "candidate_id",
                as_index=False,
                sort=False,
            )["DMS_score"]
            .mean()
        )

    frame["mutation_set"] = frame["candidate_id"].map(codec.encode)
    if frame["mutation_set"].duplicated().any():
        raise RuntimeError(
            "Reference-relative encoding produced duplicate mutation sets."
        )

    frame["crossfit_fold"] = frame["candidate_id"].map(
        crossfit_fold
    ).astype(int)
    return frame.reset_index(drop=True), codec


def deterministic_identity_subpool(
    identity_frame: pd.DataFrame,
    pool_size: int,
) -> pd.DataFrame:
    if pool_size < 1:
        raise ValueError("pool_size must be positive.")
    if len(identity_frame) < pool_size:
        raise ValueError(
            f"Assay universe has {len(identity_frame)} candidates; "
            f"requested acquisition pool {pool_size}."
        )

    work = identity_frame[
        ["candidate_id", "mutation_set", "crossfit_fold"]
    ].copy()
    work["_hash"] = work["candidate_id"].astype(str).map(
        lambda candidate_id: fnv1a32(
            f"{POOL_NAMESPACE}|{candidate_id}"
        )
    )
    return (
        work.sort_values(
            ["_hash", "candidate_id"],
            ascending=[True, True],
        )
        .head(pool_size)
        .drop(columns=["_hash"])
        .reset_index(drop=True)
    )


def bootstrap_ids(
    acquisition_pool: pd.DataFrame,
    seed_count: int,
    n_folds: int = 5,
) -> list[str]:
    if seed_count < n_folds:
        raise ValueError(
            f"seed_count must be at least {n_folds}."
        )
    if seed_count > len(acquisition_pool):
        raise ValueError("seed_count exceeds acquisition pool.")

    work = acquisition_pool[
        ["candidate_id", "crossfit_fold"]
    ].copy()
    work["_hash"] = work["candidate_id"].astype(str).map(
        lambda candidate_id: fnv1a32(
            f"{BOOTSTRAP_NAMESPACE}|{candidate_id}"
        )
    )
    work = work.sort_values(
        ["_hash", "candidate_id"],
        ascending=[True, True],
    )

    selected = []
    selected_set = set()

    for fold in range(n_folds):
        fold_rows = work[work["crossfit_fold"] == fold]
        if fold_rows.empty:
            raise RuntimeError(
                f"Acquisition pool does not populate fold {fold}."
            )
        candidate_id = str(
            fold_rows.iloc[0]["candidate_id"]
        )
        selected.append(candidate_id)
        selected_set.add(candidate_id)

    for candidate_id in work["candidate_id"].astype(str):
        if candidate_id in selected_set:
            continue
        selected.append(candidate_id)
        selected_set.add(candidate_id)
        if len(selected) == seed_count:
            break

    if len(selected) != seed_count:
        raise RuntimeError(
            "Could not build deterministic bootstrap."
        )

    return selected


def mutation_vocabulary_from_identity(
    assay_universe: pd.DataFrame,
    codec: FourSiteCodec,
) -> tuple[str, ...]:
    observed_by_position = {
        position: set()
        for position in range(1, 5)
    }

    for genotype in assay_universe["candidate_id"].astype(str):
        for position, residue in enumerate(
            genotype,
            start=1,
        ):
            observed_by_position[position].add(residue)

    tokens = []
    for position in range(1, 5):
        source = codec.reference_genotype[position - 1]
        for target in sorted(observed_by_position[position]):
            if target == source:
                continue
            tokens.append(
                f"{source}{position}{target}"
            )
    return tuple(tokens)


def sha256_ids(candidate_ids) -> str:
    payload = "\n".join(
        str(candidate_id)
        for candidate_id in candidate_ids
    )
    return hashlib.sha256(
        payload.encode("utf-8")
    ).hexdigest()
