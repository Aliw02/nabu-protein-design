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
        candidate_ids = [str(x) for x in candidate_ids]
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
        if len(set(candidate_ids)) != len(candidate_ids):
            raise ValueError("candidate_ids must be unique.")
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
        candidate_ids = [str(x) for x in candidate_ids]
        if len(mutation_sets) != len(candidate_ids):
            raise ValueError(
                "mutation_sets and candidate_ids must have equal length."
            )
        if len(set(candidate_ids)) != len(candidate_ids):
            raise ValueError("candidate_ids must be unique.")

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
            frame["V8_3_ADAPTIVE_ROUTER"] = np.nan
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
