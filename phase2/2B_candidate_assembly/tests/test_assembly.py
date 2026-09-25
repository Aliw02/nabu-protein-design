from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest


HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from assembly import (  # noqa: E402
    MeasurementAggregator,
    MutationVocabulary,
    ReferenceProtein,
    infer_reference_sequence_from_tokens,
)


def test_reference_protein_validates_source_residue():
    reference = ReferenceProtein("toy", "ACDE")
    assert reference.validate_token("A1C") == "A1C"

    with pytest.raises(ValueError):
        reference.validate_token("G1C")


def test_reference_rejects_duplicate_positions_in_candidate():
    reference = ReferenceProtein("toy", "ACDE")
    with pytest.raises(ValueError):
        reference.canonicalize_set(("A1C", "A1D"))


def test_vocabulary_is_reference_checked_and_deterministic():
    reference = ReferenceProtein("toy", "ACDE")
    vocabulary = MutationVocabulary.from_tokens(
        reference,
        ["C2A", "A1G", "A1C", "C2A"],
    )
    assert vocabulary.tokens == ("A1C", "A1G", "C2A")


def test_infer_reference_requires_consistent_sources():
    sequence = infer_reference_sequence_from_tokens(
        ["A1C", "C2A", "D3E"]
    )
    assert sequence == "ACD"

    with pytest.raises(ValueError):
        infer_reference_sequence_from_tokens(
            ["A1C", "G1D", "C2A"]
        )


def test_measurement_aggregator_freezes_mean_rule():
    reference = ReferenceProtein("toy", "ACDE")
    aggregator = MeasurementAggregator(reference)
    raw = pd.DataFrame(
        {
            "mutant": ["A1C:C2A", "C2A:A1C", "A1G"],
            "DMS_score": [1.0, 3.0, 2.0],
        }
    )

    aggregated = aggregator.aggregate(raw)
    pair = aggregated[
        aggregated["candidate_id"] == "A1C:C2A"
    ].iloc[0]

    assert pair["DMS_score"] == 2.0
    assert pair["replicate_count"] == 2
    assert np.isclose(pair["replicate_std"], np.sqrt(2.0))
