from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parents[1]
PHASE2A = HERE.parent
PHASE2A0 = PHASE2A / "2A0_campaign_simulator"
PHASE2A1 = PHASE2A / "2A1_baselines"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(PHASE2A0))
sys.path.insert(0, str(PHASE2A1))

from campaign import PolicyView, identity_pool_from_frame  # noqa: E402
from adaptive_controller import (  # noqa: E402
    EvidenceAdaptivePolicyV1,
    hierarchy_rank_disagreement,
    state_exploration_pressure,
)


def identity_pool():
    frame = pd.DataFrame(
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
    )
    return identity_pool_from_frame(frame)


def candidate_view(
    scoreable=(True, True, True),
    b3=(0.1, 0.5, 0.9),
    b4=(0.1, 0.5, 0.9),
    b5=(0.1, 0.5, 0.9),
    routed=(0.1, 0.5, 0.9),
):
    candidates = pd.DataFrame(
        {
            "candidate_id": [
                "A1C:B2D",
                "A1C:C3E",
                "B2D:C3E",
            ],
            "mutation_set": [
                ("A1C", "B2D"),
                ("A1C", "C3E"),
                ("B2D", "C3E"),
            ],
            "scoreable": scoreable,
            "B3_RAW_PAIR": b3,
            "B4_CROSSFIT_TRIPLET": b4,
            "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER": b5,
            "V8_3_ADAPTIVE_ROUTER": routed,
        }
    )
    return PolicyView(
        round_index=0,
        measurements_spent=3,
        router_mode="B3_PROTECTED_NO_HIGHER_ORDER",
        measured_ids=("A1C", "B2D", "C3E"),
        candidates=candidates,
    )


def test_pressure_is_bounded_and_monotonic_for_larger_gaps():
    low = state_exploration_pressure(0.1, 0.1, 0.1)
    high = state_exploration_pressure(0.6, 0.4, 0.3)
    assert 0.0 <= low <= 1.0
    assert 0.0 <= high <= 1.0
    assert high > low


def test_pressure_needs_no_dataset_identity():
    pressure = state_exploration_pressure(
        support_gap=0.4,
        scoreability_gap=0.2,
        disagreement_gap=0.1,
    )
    expected = 1.0 - (0.6 * 0.8 * 0.9)
    assert np.isclose(pressure, expected)


def test_hierarchy_disagreement_is_zero_when_ranks_agree():
    view = candidate_view()
    disagreement = hierarchy_rank_disagreement(view.candidates)
    assert np.allclose(disagreement.to_numpy(), 0.0)


def test_hierarchy_disagreement_increases_when_rank_order_changes():
    view = candidate_view(
        b5=(0.9, 0.1, 0.5),
    )
    disagreement = hierarchy_rank_disagreement(view.candidates)
    assert float(disagreement.max()) > 0.0


def test_controller_is_deterministic_and_records_pressure():
    pool = identity_pool()
    policy_a = EvidenceAdaptivePolicyV1(pool)
    policy_b = EvidenceAdaptivePolicyV1(pool)
    view = candidate_view()

    selected_a = policy_a.select(view, batch_size=2)
    selected_b = policy_b.select(view, batch_size=2)

    assert selected_a == selected_b
    assert policy_a.last_diagnostics is not None
    assert 0.0 <= policy_a.last_diagnostics["exploration_pressure"] <= 1.0


def test_unscoreable_candidate_receives_exploration_evidence():
    pool = identity_pool()
    policy = EvidenceAdaptivePolicyV1(pool)
    view = candidate_view(
        scoreable=(False, True, True),
        routed=(np.nan, 0.2, 0.1),
    )

    selected = policy.select(view, batch_size=1)

    assert len(selected) == 1
    assert policy.last_diagnostics["scoreability_gap"] > 0.0
