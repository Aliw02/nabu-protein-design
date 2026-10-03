from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from nabu_protein.metrics import evaluate_v83


ROOT = Path(__file__).resolve().parents[1]


def load_reference():
    path = (
        ROOT
        / "experiments"
        / "v8_3_multilandscape_validation"
        / "run_v8_3_validation.py"
    )
    spec = importlib.util.spec_from_file_location("v83_metrics_reference", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


REF = load_reference()


def test_metrics_match_frozen_reference_exactly():
    frame = pd.DataFrame(
        {
            "candidate_id": [f"C{i:03d}" for i in range(80)],
            "score": [float((i * 17) % 23) for i in range(80)],
            "truth": [float((i * 11) % 31) / 7.0 for i in range(80)],
        }
    )

    assert evaluate_v83(frame, "score", "truth") == REF.evaluate(
        frame,
        "score",
        "truth",
    )


def test_metrics_tie_breaking_is_candidate_id_deterministic():
    frame = pd.DataFrame(
        {
            "candidate_id": ["B", "A", "C"],
            "score": [1.0, 1.0, 0.0],
            "truth": [0.0, 5.0, 1.0],
        }
    )

    result = evaluate_v83(frame, "score", "truth")
    expected_pct = REF.percentile(frame["truth"].to_numpy(), 5.0)
    assert result["top1_percentile"] == expected_pct


def test_metrics_reject_nonfinite_values():
    frame = pd.DataFrame(
        {
            "candidate_id": ["A", "B"],
            "score": [1.0, np.nan],
            "truth": [1.0, 2.0],
        }
    )
    with pytest.raises(ValueError, match="finite"):
        evaluate_v83(frame, "score", "truth")
