from __future__ import annotations

from collections import Counter
from itertools import combinations
import math

import numpy as np
import pandas as pd

from nabu_protein.higher_order import fnv1a32


def _candidate_ids(frame: pd.DataFrame) -> list[str]:
    return frame["candidate_id"].astype(str).tolist()


def _score_order(
    frame: pd.DataFrame,
    score_column: str,
) -> list[str]:
    work = frame[["candidate_id", score_column]].copy()
    return (
        work.sort_values(
            [score_column, "candidate_id"],
            ascending=[False, True],
            na_position="last",
        )["candidate_id"]
        .astype(str)
        .tolist()
    )


class DeterministicRandomPolicy:
    """Fixed pseudorandom permutation, independent of assay values."""

    name = "deterministic_random"

    def __init__(self, seed: int = 161):
        self.seed = int(seed)

    def select(self, view, batch_size: int) -> list[str]:
        work = view.candidates[["candidate_id"]].copy()
        work["_random_key"] = work["candidate_id"].map(
            lambda candidate_id: fnv1a32(
                f"NABU_PHASE2A1_RANDOM|{self.seed}|{candidate_id}"
            )
        )
        work = work.sort_values(
            ["_random_key", "candidate_id"],
            ascending=[True, True],
        )
        return work.head(batch_size)["candidate_id"].astype(str).tolist()


class GreedyV83Policy:
    """Exploit the current V8.3 within-pool ordering."""

    name = "greedy_v83"

    def select(self, view, batch_size: int) -> list[str]:
        order = _score_order(
            view.candidates,
            "V8_3_ADAPTIVE_ROUTER",
        )
        return order[:batch_size]


class StaticInitialV83Policy:
    """Freeze the V8.3 ranking from the first post-bootstrap fit."""

    name = "static_initial_v83"

    def __init__(self):
        self._frozen_order: list[str] | None = None

    def select(self, view, batch_size: int) -> list[str]:
        if self._frozen_order is None:
            self._frozen_order = _score_order(
                view.candidates,
                "V8_3_ADAPTIVE_ROUTER",
            )

        available = set(_candidate_ids(view.candidates))
        selected = [
            candidate_id
            for candidate_id in self._frozen_order
            if candidate_id in available
        ][:batch_size]

        if len(selected) != batch_size:
            raise RuntimeError(
                "Static initial ranking could not supply the requested batch."
            )
        return selected


class _SupportPolicyBase:
    def __init__(self, identity_pool: pd.DataFrame):
        lookup = identity_pool[
            ["candidate_id", "mutation_set"]
        ].copy()
        lookup["candidate_id"] = lookup["candidate_id"].map(str)
        self._mutation_by_id = dict(
            zip(
                lookup["candidate_id"].tolist(),
                lookup["mutation_set"].tolist(),
            )
        )

    def structural_exploration(self, view) -> pd.Series:
        main_support: Counter[str] = Counter()
        pair_support: Counter[tuple[str, str]] = Counter()

        for candidate_id in view.measured_ids:
            mutation_set = self._mutation_by_id[str(candidate_id)]
            main_support.update(mutation_set)
            pair_support.update(combinations(mutation_set, 2))

        values = []
        for mutation_set in view.candidates["mutation_set"]:
            supports = [
                int(main_support[mutation])
                for mutation in mutation_set
            ]
            supports.extend(
                int(pair_support[pair])
                for pair in combinations(mutation_set, 2)
            )
            if not supports:
                values.append(1.0)
                continue
            deficits = [
                1.0 / math.sqrt(1.0 + float(support))
                for support in supports
            ]
            values.append(float(np.mean(deficits)))

        return pd.Series(values, index=view.candidates.index, dtype=float)


class ExplorationOnlyPolicy(_SupportPolicyBase):
    """Historical main/pair support-deficit exploration only."""

    name = "exploration_only"

    def select(self, view, batch_size: int) -> list[str]:
        work = view.candidates[["candidate_id"]].copy()
        work["structural_exploration"] = self.structural_exploration(view)
        work = work.sort_values(
            ["structural_exploration", "candidate_id"],
            ascending=[False, True],
        )
        return work.head(batch_size)["candidate_id"].astype(str).tolist()


class HistoricalFiftyFiftyPolicy(_SupportPolicyBase):
    """Frozen historical 0.5 exploration + 0.5 exploitation rank fusion."""

    name = "historical_50_50"

    def select(self, view, batch_size: int) -> list[str]:
        work = view.candidates[
            ["candidate_id", "V8_3_ADAPTIVE_ROUTER"]
        ].copy()

        work["structural_exploration"] = self.structural_exploration(view)
        work["exploration_rank"] = work[
            "structural_exploration"
        ].rank(method="average", pct=True)

        work["exploitation_rank"] = work[
            "V8_3_ADAPTIVE_ROUTER"
        ].rank(method="average", pct=True).fillna(0.0)

        work["acquisition_score"] = (
            0.5 * work["exploration_rank"]
            + 0.5 * work["exploitation_rank"]
        )

        work = work.sort_values(
            ["acquisition_score", "candidate_id"],
            ascending=[False, True],
        )
        return work.head(batch_size)["candidate_id"].astype(str).tolist()


def build_baselines(identity_pool: pd.DataFrame):
    return [
        DeterministicRandomPolicy(seed=161),
        GreedyV83Policy(),
        ExplorationOnlyPolicy(identity_pool),
        HistoricalFiftyFiftyPolicy(identity_pool),
        StaticInitialV83Policy(),
    ]
