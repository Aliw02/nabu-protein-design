from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd


HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from run_multilandscape_evaluation import (  # noqa: E402
    POOL_NAMESPACE,
    deterministic_identity_subpool,
    load_standard_landscape,
    parse_landscape,
)


def test_identity_subpool_is_deterministic_and_truth_independent():
    frame = pd.DataFrame(
        {
            "candidate_id": [f"X1A:X2A:X3A:X4{i}" for i in range(20)],
            "mutant": [f"X1A:X2A:X3A:X4{i}" for i in range(20)],
            "DMS_score": [float(i) for i in range(20)],
        }
    )
    altered = frame.copy()
    altered["DMS_score"] = list(reversed(altered["DMS_score"].tolist()))

    left = deterministic_identity_subpool(frame, pool_size=10)
    right = deterministic_identity_subpool(altered, pool_size=10)

    assert left["candidate_id"].tolist() == right["candidate_id"].tolist()
    assert len(left) == 10
    assert POOL_NAMESPACE


def test_standard_landscape_loader_creates_candidate_ids(tmp_path):
    path = tmp_path / "landscape.csv"
    pd.DataFrame(
        {
            "mutant": ["X1A:X2C:X3D:X4E", "X1C:X2D:X3E:X4F"],
            "DMS_score": [0.1, 0.2],
        }
    ).to_csv(path, index=False)

    loaded = load_standard_landscape(path)

    assert loaded["candidate_id"].tolist() == loaded["mutant"].tolist()
    assert len(loaded) == 2


def test_duplicate_candidates_are_aggregated_before_campaign():
    frame = pd.DataFrame(
        {
            "mutant": ["X1A:X2C:X3D:X4E", "X1A:X2C:X3D:X4E"],
            "DMS_score": [0.1, 0.3],
        }
    )
    path = HERE / "_tmp_duplicate_landscape.csv"
    try:
        frame.to_csv(path, index=False)
        loaded = load_standard_landscape(path)
        assert len(loaded) == 1
        assert loaded.iloc[0]["DMS_score"] == 0.2
    finally:
        if path.exists():
            path.unlink()


def test_landscape_cli_parser():
    assert parse_landscape("GB1=/tmp/gb1.csv") == (
        "GB1",
        "/tmp/gb1.csv",
    )
