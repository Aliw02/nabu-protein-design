from __future__ import annotations

import gzip
from pathlib import Path
import sys

import pandas as pd

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from ired_scoreability_preflight import scoreability_preflight


def _write(path: Path, rows: list[dict]) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        pd.DataFrame(rows).to_csv(handle, index=False)


def test_scoreability_uses_identity_only(tmp_path):
    path = tmp_path / "ired.csv.gz"
    _write(
        path,
        [
            {"sequence": "AAAA", "target": 999.0, "set": "train", "validation": False},
            {"sequence": "CAAA", "target": -999.0, "set": "train", "validation": False},
            {"sequence": "ACAA", "target": 123.0, "set": "train", "validation": False},
            {"sequence": "CCAA", "target": 456.0, "set": "train", "validation": False},
            {"sequence": "CCCA", "target": -1.0, "set": "test", "validation": False},
        ],
    )
    result = scoreability_preflight(path)
    assert result["target_values_read"] is False
    assert result["fit_count"] == 4
    assert result["test_count"] == 1


def test_unseen_main_makes_test_unscoreable(tmp_path):
    path = tmp_path / "ired.csv.gz"
    _write(
        path,
        [
            {"sequence": "AAAA", "target": 0.0, "set": "train", "validation": False},
            {"sequence": "CAAA", "target": 0.0, "set": "train", "validation": False},
            {"sequence": "ACAA", "target": 0.0, "set": "train", "validation": False},
            {"sequence": "CCAA", "target": 0.0, "set": "train", "validation": False},
            {"sequence": "CCDA", "target": 1.0, "set": "test", "validation": False},
        ],
    )
    result = scoreability_preflight(path)
    assert result["scoreable_count"] == 0
    assert result["unscoreable_count"] == 1
    assert not result["full_test_scoreable"]
