"""
Reusable facade for the frozen NABU V8.3 Phase-1 core.
"""

from __future__ import annotations

from itertools import combinations
import re

import numpy as np
import pandas as pd

from .higher_order import fit_crossfitted_hierarchy, score_hierarchy
from .router import apply_router, decide_router


_MUTATION_TOKEN_RE = re.compile(r"[A-Z][0-9]+[A-Z*]")

_SCORE_COLUMNS = (
    "B2_ADDITIVE",
    "B3_RAW_PAIR",
    "B4_CROSSFIT_TRIPLET",
    "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER",
    "triplet_delta",
    "quartet_delta",
    "triplet_supported",
    "quartet_supported",
    "triplet_confidence",
    "quartet_confidence",
    "V8_3_ADAPTIVE_ROUTER",
)


def _normalize_candidate_ids(candidate_ids):
    normalized = []
    for value in candidate_ids:
        if value is None:
            raise ValueError("candidate_ids must not contain None values.")
        # Preserve the historical Phase-1 behavior for CSV blank identities:
        # pandas NaN stringifies to "nan" before deterministic fold hashing.
        text = str(value)
        if not text.strip():
            raise ValueError("candidate_ids must not contain blank values.")
        normalized.append(text)
    if len(set(normalized)) != len(normalized):
        raise ValueError("candidate_ids must be unique.")
    return normalized


def canonicalize_mutation_set(mutation_set):
    """Validate and canonicalize one pre-parsed mutation collection."""
    if isinstance(mutation_set, str):
        raise ValueError(
            "mutation_set must be an iterable of mutation tokens, not a raw string."
        )
    tokens = tuple(str(token) for token in mutation_set)
    invalid = [
        token
        for token in tokens
        if _MUTATION_TOKEN_RE.fullmatch(token) is None
    ]
    if invalid:
        raise ValueError(f"Invalid mutation token(s): {invalid}")
    if len(set(tokens)) != len(tokens):
        raise ValueError("mutation_set contains duplicate mutation tokens.")

    parsed = [
        (int(token[1:-1]), token[0], token[-1], token)
        for token in tokens
    ]
    if any(position < 1 for position, _, _, _ in parsed):
        raise ValueError("mutation positions must be 1-based positive integers.")

    positions = [position for position, _, _, _ in parsed]
    if len(set(positions)) != len(positions):
        raise ValueError(
            "mutation_set contains multiple substitutions at the same residue position."
        )

    no_ops = [
        token
        for _, source, target, token in parsed
        if target != "*" and source == target
    ]
    if no_ops:
        raise ValueError(f"Mutation token(s) do not change residue state: {no_ops}")

    return tuple(
        sorted(
            tokens,
            key=lambda token: (int(token[1:-1]), token),
        )
    )


def is_scoreable(mutation_set, model) -> bool:
    mutation_set = canonicalize_mutation_set(mutation_set)
    return bool(
        all(
            mutation in model["base"]["main"]
            for mutation in mutation_set
        )
        and any(
            pair in model["base"]["pair"]
            for pair in combinations(mutation_set, 2)
        )
    )


class NabuV83Model:
    """Frozen B2/B3/B4/B5 hierarchy plus the V8.3 OOF router."""

    def __init__(self):
        self.model = None
        self.visible_ids = None
        self.visible_labels = None
        self.router_decision = None

    def fit(self, mutation_sets, labels, candidate_ids):
        mutation_sets = list(mutation_sets)
        candidate_ids = _normalize_candidate_ids(candidate_ids)
        labels = np.asarray(labels, dtype=float)

        if not (
            len(mutation_sets) == len(labels) == len(candidate_ids)
        ):
            raise ValueError(
                "mutation_sets, labels, and candidate_ids must have equal length."
            )
        if len(labels) == 0:
            raise ValueError("At least one visible candidate is required.")
        if not np.isfinite(labels).all():
            raise ValueError("labels must contain only finite numeric values.")
        canonical_sets = [
            canonicalize_mutation_set(ms)
            for ms in mutation_sets
        ]
        if len(set(canonical_sets)) != len(canonical_sets):
            raise ValueError(
                "mutation_sets must be unique; aggregate or deduplicate replicate rows before fit."
            )

        self.model = fit_crossfitted_hierarchy(
            canonical_sets,
            labels,
            candidate_ids,
        )
        self.visible_ids = candidate_ids
        self.visible_labels = labels

        self.router_decision = decide_router(
            self.model,
            self.visible_ids,
            self.visible_labels,
        )
        return self

    def score_candidates(self, mutation_sets, candidate_ids):
        if self.model is None:
            raise RuntimeError("Model must be fit before scoring.")

        mutation_sets = [
            canonicalize_mutation_set(ms)
            for ms in mutation_sets
        ]
        candidate_ids = _normalize_candidate_ids(candidate_ids)
        if len(mutation_sets) != len(candidate_ids):
            raise ValueError(
                "mutation_sets and candidate_ids must have equal length."
            )
        frame = pd.DataFrame(
            {
                "candidate_id": candidate_ids,
                "mutation_set": mutation_sets,
            }
        )
        frame["scoreable"] = frame["mutation_set"].map(
            lambda mutation_set: is_scoreable(
                mutation_set,
                self.model,
            )
        )

        scoreable = frame[frame["scoreable"]].copy()
        if scoreable.empty:
            for column in _SCORE_COLUMNS:
                frame[column] = np.nan
            return frame

        rows = [
            score_hierarchy(mutation_set, self.model)
            for mutation_set in scoreable["mutation_set"]
        ]
        scores = pd.DataFrame(rows, index=scoreable.index)
        for column in scores.columns:
            scoreable[column] = scores[column]

        apply_router(
            scoreable,
            self.model,
            self.visible_ids,
            self.visible_labels,
            "B3_RAW_PAIR",
            "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER",
            "V8_3_ADAPTIVE_ROUTER",
        )

        score_columns = list(scores.columns) + [
            "V8_3_ADAPTIVE_ROUTER"
        ]
        for column in score_columns:
            frame[column] = np.nan
            frame.loc[scoreable.index, column] = scoreable[column]

        return frame
