"""
Reusable facade for the frozen NABU V8.3 Phase-1 core.
"""

from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd

from .higher_order import fit_crossfitted_hierarchy, score_hierarchy
from .router import apply_router


def is_scoreable(mutation_set, model) -> bool:
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
        if len(set(candidate_ids)) != len(candidate_ids):
            raise ValueError("candidate_ids must be unique.")

        self.model = fit_crossfitted_hierarchy(
            mutation_sets,
            labels,
            candidate_ids,
        )
        self.visible_ids = candidate_ids
        self.visible_labels = labels

        dummy = pd.DataFrame(
            {
                "candidate_id": ["DUMMY"],
                "B3_RAW_PAIR": [0.0],
                "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER": [0.0],
            }
        )
        self.router_decision = apply_router(
            dummy,
            self.model,
            self.visible_ids,
            self.visible_labels,
            "B3_RAW_PAIR",
            "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER",
            "V8_3_ADAPTIVE_ROUTER",
        )
        return self

    def score_candidates(self, mutation_sets, candidate_ids):
        if self.model is None:
            raise RuntimeError("Model must be fit before scoring.")

        mutation_sets = list(mutation_sets)
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
