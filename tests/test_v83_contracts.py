from __future__ import annotations

import pandas as pd

from nabu_protein.higher_order import parse_mutations
from nabu_protein.router import rank_preserving_elite_rerank
from nabu_protein.v83 import is_scoreable


def test_canonical_parser_orders_by_position_and_keeps_stop_token():
    assert parse_mutations("D4F:A1C:C3*") == (
        "A1C",
        "C3*",
        "D4F",
    )


def test_rank_preserving_router_never_moves_outside_candidate_into_elite():
    frame = pd.DataFrame(
        {
            "candidate_id": ["B", "A", "C", "D"],
            "B3": [10.0, 10.0, 9.0, 8.0],
            "B5": [2.0, 1.0, 100.0, 90.0],
        }
    )

    details = rank_preserving_elite_rerank(
        frame,
        "B3",
        "B5",
        "OUT",
        fraction=0.50,
    )

    ranked = frame.sort_values(
        ["OUT", "candidate_id"],
        ascending=[False, True],
    )["candidate_id"].tolist()

    assert ranked == ["B", "A", "C", "D"]
    assert details["elite_count"] == 2
    assert details["boundary_preserved"] is True
    assert details["outside_B3_order_preserved"] is True


def test_scoreability_matches_frozen_phase1_rule():
    model = {
        "base": {
            "main": {
                "A1C": {},
                "B2D": {},
                "C3E": {},
            },
            "pair": {
                ("A1C", "B2D"): {},
            },
        }
    }

    assert is_scoreable(("A1C", "B2D"), model) is True
    assert is_scoreable(("A1C", "B2D", "C3E"), model) is True
    assert is_scoreable(("A1C", "C3E"), model) is False
    assert is_scoreable(("A1C", "D4F"), model) is False


def test_facade_rejects_duplicate_candidates_and_nonfinite_labels():
    import numpy as np
    import pytest

    from nabu_protein.v83 import NabuV83Model

    model = NabuV83Model()

    with pytest.raises(ValueError, match="candidate_ids must be unique"):
        model.fit(
            [("A1C", "B2D"), ("A1C", "C3E")],
            [1.0, 2.0],
            ["dup", "dup"],
        )

    with pytest.raises(ValueError, match="labels must contain only finite"):
        model.fit(
            [("A1C", "B2D"), ("A1C", "C3E")],
            [1.0, np.nan],
            ["a", "b"],
        )


def test_facade_rejects_duplicate_mutation_rows_before_crossfit():
    import pytest

    from nabu_protein.v83 import NabuV83Model

    with pytest.raises(ValueError, match="mutation_sets must be unique"):
        NabuV83Model().fit(
            [("A1C", "B2D"), ("A1C", "B2D")],
            [1.0, 1.1],
            ["replicate-a", "replicate-b"],
        )
