from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parents[1]
PHASE2A0 = HERE.parent / "2A0_campaign_simulator"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(PHASE2A0))

from campaign import PolicyView, identity_pool_from_frame  # noqa: E402
from baselines import (  # noqa: E402
    DeterministicRandomPolicy,
    ExplorationOnlyPolicy,
    GreedyV83Policy,
    HistoricalFiftyFiftyPolicy,
    StaticInitialV83Policy,
)


def view(scores, measured_ids=()):
    candidates = pd.DataFrame(
        {
            "candidate_id": ["A1C:B2D", "A1C:C3E", "B2D:C3E"],
            "mutation_set": [
                ("A1C", "B2D"),
                ("A1C", "C3E"),
                ("B2D", "C3E"),
            ],
            "V8_3_ADAPTIVE_ROUTER": scores,
            "scoreable": [True, True, True],
        }
    )
    return PolicyView(
        round_index=0,
        measurements_spent=len(measured_ids),
        router_mode="B3_PROTECTED_NO_HIGHER_ORDER",
        measured_ids=tuple(measured_ids),
        candidates=candidates,
    )


def identity_pool():
    return pd.DataFrame(
        {
            "candidate_id": [
                "A1C",
                "B2D",
                "C3E",
                "A1C:B2D",
                "A1C:C3E",
                "B2D:C3E",
            ],
            "mutant": [
                "A1C",
                "B2D",
                "C3E",
                "A1C:B2D",
                "A1C:C3E",
                "B2D:C3E",
            ],
        }
    ).pipe(identity_pool_from_frame)


def test_greedy_uses_current_v83_order():
    policy = GreedyV83Policy()
    selected = policy.select(view([0.2, 0.9, 0.5]), batch_size=2)
    assert selected == ["A1C:C3E", "B2D:C3E"]


def test_static_initial_ranking_does_not_follow_later_score_changes():
    policy = StaticInitialV83Policy()
    first = policy.select(view([0.9, 0.5, 0.2]), batch_size=1)
    assert first == ["A1C:B2D"]

    second_view = view([0.1, 0.2, 9.0])
    second_view = PolicyView(
        round_index=1,
        measurements_spent=1,
        router_mode=second_view.router_mode,
        measured_ids=("A1C:B2D",),
        candidates=second_view.candidates[
            second_view.candidates["candidate_id"] != "A1C:B2D"
        ].copy(),
    )
    second = policy.select(second_view, batch_size=1)
    assert second == ["A1C:C3E"]


def test_deterministic_random_replays_identically():
    policy_a = DeterministicRandomPolicy(seed=161)
    policy_b = DeterministicRandomPolicy(seed=161)
    current = view([0.1, 0.2, 0.3])
    assert policy_a.select(current, 3) == policy_b.select(current, 3)


def test_exploration_prefers_lower_supported_structure():
    pool = identity_pool()
    policy = ExplorationOnlyPolicy(pool)

    current = view(
        [0.0, 0.0, 0.0],
        measured_ids=("A1C", "B2D", "A1C:B2D"),
    )
    selected = policy.select(current, 1)

    assert selected[0] in {"A1C:C3E", "B2D:C3E"}


def test_historical_50_50_returns_exact_batch_without_truth():
    pool = identity_pool()
    policy = HistoricalFiftyFiftyPolicy(pool)
    current = view(
        [0.2, np.nan, 0.8],
        measured_ids=("A1C", "B2D", "C3E"),
    )

    selected = policy.select(current, 2)

    assert len(selected) == 2
    assert len(set(selected)) == 2
    assert set(selected).issubset(set(current.candidates["candidate_id"]))
