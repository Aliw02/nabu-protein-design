from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
import re

import numpy as np
import pandas as pd

from nabu_protein.v83 import canonicalize_mutation_set


STANDARD_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
TOKEN_RE = re.compile(r"([A-Z])([1-9][0-9]*)([A-Z])")


@dataclass(frozen=True)
class ReferenceProtein:
    name: str
    sequence: str

    def __post_init__(self):
        sequence = str(self.sequence).strip().upper()
        if not sequence:
            raise ValueError("Reference sequence must not be blank.")
        invalid = sorted(set(sequence) - STANDARD_AA)
        if invalid:
            raise ValueError(
                f"Reference sequence contains non-standard residues: {invalid}"
            )
        object.__setattr__(self, "sequence", sequence)

    def residue_at(self, position: int) -> str:
        position = int(position)
        if position < 1 or position > len(self.sequence):
            raise ValueError(
                f"Position {position} outside reference length "
                f"{len(self.sequence)}."
            )
        return self.sequence[position - 1]

    def validate_token(self, token: str) -> str:
        token = str(token).strip().upper()
        match = TOKEN_RE.fullmatch(token)
        if match is None:
            raise ValueError(f"Invalid mutation token: {token!r}")

        source, position_text, target = match.groups()
        position = int(position_text)

        expected = self.residue_at(position)
        if source != expected:
            raise ValueError(
                f"Mutation {token} source residue {source} does not match "
                f"reference {expected} at position {position}."
            )
        if target not in STANDARD_AA:
            raise ValueError(
                f"Mutation {token} target is not a standard amino acid."
            )
        if target == source:
            raise ValueError(f"Mutation {token} is a no-op.")
        return token

    def canonicalize_set(self, mutation_set) -> tuple[str, ...]:
        tokens = tuple(
            self.validate_token(token)
            for token in mutation_set
        )
        canonical = canonicalize_mutation_set(tokens)
        return canonical


@dataclass(frozen=True)
class MutationVocabulary:
    reference: ReferenceProtein
    tokens: tuple[str, ...]

    @classmethod
    def from_tokens(
        cls,
        reference: ReferenceProtein,
        tokens,
    ) -> "MutationVocabulary":
        canonical_tokens = sorted(
            {
                reference.validate_token(token)
                for token in tokens
            },
            key=lambda token: (int(token[1:-1]), token),
        )
        if not canonical_tokens:
            raise ValueError("Mutation vocabulary must not be empty.")
        return cls(reference=reference, tokens=tuple(canonical_tokens))

    @property
    def positions(self) -> tuple[int, ...]:
        return tuple(
            sorted({int(token[1:-1]) for token in self.tokens})
        )


def infer_reference_sequence_from_tokens(tokens) -> str:
    source_by_position: dict[int, str] = {}

    for raw in tokens:
        token = str(raw).strip().upper()
        match = TOKEN_RE.fullmatch(token)
        if match is None:
            raise ValueError(f"Invalid mutation token: {token!r}")
        source, position_text, _ = match.groups()
        position = int(position_text)

        previous = source_by_position.get(position)
        if previous is not None and previous != source:
            raise ValueError(
                f"Inconsistent source residues at position {position}: "
                f"{previous} vs {source}."
            )
        source_by_position[position] = source

    if not source_by_position:
        raise ValueError("Cannot infer reference from an empty vocabulary.")

    max_position = max(source_by_position)
    missing = [
        position
        for position in range(1, max_position + 1)
        if position not in source_by_position
    ]
    if missing:
        raise ValueError(
            "Development reference inference requires contiguous positions; "
            f"missing {missing}."
        )

    return "".join(
        source_by_position[position]
        for position in range(1, max_position + 1)
    )


class MeasurementAggregator:
    def __init__(self, reference: ReferenceProtein):
        self.reference = reference

    def aggregate(
        self,
        frame: pd.DataFrame,
        mutation_column: str = "mutant",
        label_column: str = "DMS_score",
    ) -> pd.DataFrame:
        required = {mutation_column, label_column}
        missing = required - set(frame.columns)
        if missing:
            raise ValueError(
                f"Measurement frame missing columns: {sorted(missing)}"
            )

        rows = []
        for row in frame[[mutation_column, label_column]].itertuples(index=False):
            raw_mutant = str(row[0]).strip()
            if not raw_mutant:
                raise ValueError("Measurement mutation identity must not be blank.")
            tokens = tuple(
                token.strip()
                for token in raw_mutant.split(":")
                if token.strip()
            )
            canonical = self.reference.canonicalize_set(tokens)
            label = float(row[1])
            if not np.isfinite(label):
                raise ValueError("Measurement labels must be finite.")
            rows.append(
                {
                    "mutation_set": canonical,
                    "candidate_id": ":".join(canonical),
                    "label": label,
                }
            )

        work = pd.DataFrame(rows)
        grouped = (
            work.groupby(
                ["candidate_id", "mutation_set"],
                as_index=False,
                sort=False,
            )["label"]
            .agg(["mean", "count", "std"])
            .reset_index()
        )
        grouped = grouped.rename(
            columns={
                "mean": "DMS_score",
                "count": "replicate_count",
                "std": "replicate_std",
            }
        )
        grouped["replicate_std"] = grouped[
            "replicate_std"
        ].fillna(0.0)
        grouped["mutant"] = grouped["candidate_id"]
        return grouped[
            [
                "candidate_id",
                "mutant",
                "mutation_set",
                "DMS_score",
                "replicate_count",
                "replicate_std",
            ]
        ]


def _candidate_id(mutation_set: tuple[str, ...]) -> str:
    return ":".join(mutation_set)


class CandidateAssembler:
    def __init__(
        self,
        model,
        vocabulary: MutationVocabulary,
        measured_mutation_sets,
    ):
        if model.model is None:
            raise ValueError("NABU V8.3 model must be fitted.")

        self.model = model
        self.vocabulary = vocabulary
        self.measured = {
            vocabulary.reference.canonicalize_set(mutation_set)
            for mutation_set in measured_mutation_sets
        }

    def _valid_expand(
        self,
        partial: tuple[str, ...],
        token: str,
    ) -> tuple[str, ...] | None:
        position = int(token[1:-1])
        used_positions = {
            int(existing[1:-1])
            for existing in partial
        }
        if position in used_positions:
            return None

        candidate = self.vocabulary.reference.canonicalize_set(
            tuple(partial) + (token,)
        )
        if candidate in self.measured:
            return None
        return candidate

    def _score(
        self,
        candidates: list[tuple[str, ...]],
    ) -> pd.DataFrame:
        candidate_ids = [_candidate_id(candidate) for candidate in candidates]
        scored = self.model.score_candidates(
            candidates,
            candidate_ids,
        )
        scored["mutation_count"] = scored[
            "mutation_set"
        ].map(len)
        return scored

    @staticmethod
    def _rank(scored: pd.DataFrame) -> pd.DataFrame:
        work = scored.copy()
        work["_scoreable_sort"] = work["scoreable"].astype(int)
        work["_router_sort"] = work[
            "V8_3_ADAPTIVE_ROUTER"
        ].astype(float).fillna(-np.inf)
        return work.sort_values(
            ["_scoreable_sort", "_router_sort", "candidate_id"],
            ascending=[False, False, True],
        ).drop(columns=["_scoreable_sort", "_router_sort"])

    def assemble(
        self,
        target_order: int = 3,
        beam_width: int = 256,
        proposal_count: int = 64,
    ) -> tuple[pd.DataFrame, list[dict]]:
        if target_order < 2:
            raise ValueError("target_order must be at least 2.")
        if beam_width < 1 or proposal_count < 1:
            raise ValueError("beam_width and proposal_count must be positive.")

        pair_candidates = []
        for left, right in combinations(self.vocabulary.tokens, 2):
            candidate = self._valid_expand((left,), right)
            if candidate is not None:
                pair_candidates.append(candidate)

        pair_candidates = sorted(set(pair_candidates))
        pair_scored = self._score(pair_candidates)
        pair_ranked = self._rank(pair_scored)

        trace = [
            {
                "order": 2,
                "generated_unique": int(len(pair_candidates)),
                "scoreable": int(pair_scored["scoreable"].sum()),
                "beam_kept": int(min(beam_width, len(pair_ranked))),
            }
        ]

        if target_order == 2:
            proposals = pair_ranked.head(proposal_count).copy()
            return proposals.reset_index(drop=True), trace

        frontier = [
            tuple(value)
            for value in pair_ranked.head(beam_width)["mutation_set"]
        ]

        for order in range(3, target_order + 1):
            expanded = set()
            for partial in frontier:
                for token in self.vocabulary.tokens:
                    candidate = self._valid_expand(partial, token)
                    if candidate is None:
                        continue
                    if len(candidate) != order:
                        continue
                    expanded.add(candidate)

            candidates = sorted(expanded)
            if not candidates:
                raise RuntimeError(
                    f"Assembly frontier became empty at order {order}."
                )

            scored = self._score(candidates)
            ranked = self._rank(scored)

            trace.append(
                {
                    "order": int(order),
                    "generated_unique": int(len(candidates)),
                    "scoreable": int(scored["scoreable"].sum()),
                    "beam_kept": int(min(beam_width, len(ranked))),
                }
            )

            if order == target_order:
                proposals = ranked.head(proposal_count).copy()
                return proposals.reset_index(drop=True), trace

            frontier = [
                tuple(value)
                for value in ranked.head(beam_width)["mutation_set"]
            ]

        raise RuntimeError("Assembly did not produce a proposal set.")
