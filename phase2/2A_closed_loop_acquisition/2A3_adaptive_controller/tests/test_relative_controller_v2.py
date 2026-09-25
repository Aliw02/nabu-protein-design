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
    RelativeEvidencePolicyV2,
    relative_exploration_pressure,
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
    routed=(0.1, 0.5, 0.9),
    scoreable=(True, True, True),
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
            "B3_RAW_PAIR": routed,
            "B4_CROSSFIT_TRIPLET": routed,
            "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER": routed,
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


def test_relative_pressure_is_neutral_at_reference():
    assert np.isclose(
        relative_exploration_pressure(0.6, 0.6),
        0.5,
    )


def test_relative_pressure_falls_when_gap_improves():
    improved = relative_exploration_pressure(0.3, 0.6)
    neutral = relative_exploration_pressure(0.6, 0.6)
    assert improved < neutral


def test_relative_pressure_rises_when_gap_worsens():
    worsened = relative_exploration_pressure(0.8, 0.6)
    neutral = relative_exploration_pressure(0.6, 0.6)
    assert worsened > neutral


def test_zero_gap_has_zero_pressure():
    assert relative_exploration_pressure(0.0, 0.0) == 0.0


def test_v2_initial_round_sets_reference_and_starts_at_half():
    policy = RelativeEvidencePolicyV2(identity_pool())
    selected = policy.select(candidate_view(), batch_size=2)

    assert len(selected) == 2
    assert policy.reference_gap is not None
    assert np.isclose(
        policy.last_diagnostics["exploration_pressure"],
        0.5,
    )


def test_v2_is_deterministic():
    pool = identity_pool()
    view = candidate_view()
    a = RelativeEvidencePolicyV2(pool)
    b = RelativeEvidencePolicyV2(pool)

    assert a.select(view, 2) == b.select(view, 2)
