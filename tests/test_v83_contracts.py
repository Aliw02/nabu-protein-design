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


def test_public_mutation_sets_are_canonicalized_before_lookup():
    from nabu_protein.v83 import canonicalize_mutation_set

    assert canonicalize_mutation_set(("C3E", "A1C", "B2D")) == (
        "A1C",
        "B2D",
        "C3E",
    )


def test_public_mutation_set_validation_rejects_malformed_and_duplicates():
    import pytest

    from nabu_protein.v83 import canonicalize_mutation_set

    with pytest.raises(ValueError, match="raw string"):
        canonicalize_mutation_set("A1C:B2D")

    with pytest.raises(ValueError, match="Invalid mutation"):
        canonicalize_mutation_set(("A1C", "bad"))

    with pytest.raises(ValueError, match="duplicate mutation"):
        canonicalize_mutation_set(("A1C", "A1C"))


def test_public_mutation_set_validation_rejects_impossible_same_position_and_noop():
    import pytest

    from nabu_protein.v83 import canonicalize_mutation_set

    with pytest.raises(ValueError, match="same residue position"):
        canonicalize_mutation_set(("A10C", "A10D"))

    with pytest.raises(ValueError, match="do not change residue state"):
        canonicalize_mutation_set(("A10A",))

    with pytest.raises(ValueError, match="1-based positive"):
        canonicalize_mutation_set(("A0C",))


def test_facade_rejects_null_and_blank_candidate_ids():
    import pytest

    from nabu_protein.v83 import NabuV83Model

    with pytest.raises(ValueError, match="null values"):
        NabuV83Model().fit(
            [("A1C", "B2D"), ("A1C", "C3E")],
            [1.0, 2.0],
            ["a", None],
        )

    with pytest.raises(ValueError, match="blank values"):
        NabuV83Model().fit(
            [("A1C", "B2D"), ("A1C", "C3E")],
            [1.0, 2.0],
            ["a", "   "],
        )


def test_score_candidates_empty_scoreable_set_keeps_stable_schema():
    from nabu_protein.v83 import NabuV83Model

    candidate_ids = [
        "V0",
        "V1",
        "V2",
        "V3",
        "V4",
        "V5",
        "V6",
        "V7",
        "V8",
        "V9",
    ]
    mutation_sets = [
        ("A1C", "B2D"),
        ("A1C", "C3E"),
        ("A1C", "D4F"),
        ("B2D", "C3E"),
        ("B2D", "D4F"),
        ("C3E", "D4F"),
        ("A1C", "B2D", "C3E"),
        ("A1C", "B2D", "D4F"),
        ("A1C", "C3E", "D4F"),
        ("B2D", "C3E", "D4F"),
    ]
    labels = [float(i) for i in range(len(candidate_ids))]

    model = NabuV83Model().fit(
        mutation_sets,
        labels,
        candidate_ids,
    )
    scored = model.score_candidates(
        [("X20Y", "Z21A")],
        ["unseen"],
    )

    expected_columns = {
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
    }
    assert expected_columns.issubset(scored.columns)
    assert scored["scoreable"].tolist() == [False]
    assert scored[list(expected_columns)].isna().all().all()
