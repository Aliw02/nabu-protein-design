from __future__ import annotations

import gzip
from pathlib import Path
import sys
import zipfile

import pandas as pd

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from preflight_identity import preflight


def _write_fixture(tmp_path: Path):
    aav_zip = tmp_path / "aav.zip"
    aav_csv = pd.DataFrame(
        {
            "sequence": ["ACDE", "ACDF"],
            "target": [999.0, -999.0],
            "set": ["train", "test"],
            "validation": [None, False],
        }
    )
    with zipfile.ZipFile(aav_zip, "w") as archive:
        archive.writestr(
            "splits/regression/sampled.csv",
            aav_csv.to_csv(index=False),
        )

    fasta = tmp_path / "wt.fasta"
    fasta.write_text(">wt\nACDE\n", encoding="utf-8")

    ired = tmp_path / "ired.csv.gz"
    ired_frame = pd.DataFrame(
        {
            "sequence": ["ACDE", "ACDF", "ACDG"],
            "target": [1000.0, -1000.0, 0.0],
            "set": ["train", "train", "test"],
            "validation": [False, True, False],
        }
    )
    with gzip.open(ired, "wt", encoding="utf-8") as handle:
        ired_frame.to_csv(handle, index=False)

    return aav_zip, fasta, ired


def test_preflight_is_identity_only_and_hashes_sources(tmp_path):
    aav_zip, fasta, ired = _write_fixture(tmp_path)
    result = preflight(aav_zip, fasta, ired)
    assert result["target_values_read"] is False
    assert len(result["sources"]["aav_splits_zip_sha256"]) == 64
    assert len(result["sources"]["ired_two_to_many_gz_sha256"]) == 64


def test_aav_gate_requires_complete_representability(tmp_path):
    aav_zip, fasta, ired = _write_fixture(tmp_path)
    result = preflight(aav_zip, fasta, ired)
    assert result["gates"]["aav_identity_gate_pass"]


def test_ired_reference_is_derived_without_target(tmp_path):
    aav_zip, fasta, ired = _write_fixture(tmp_path)
    result = preflight(aav_zip, fasta, ired)
    assert result["ired"]["derived_reference"] == "ACDE"
    assert result["ired"]["reference_present_as_exact_row"]
    assert result["ired"]["fit_train_rows"] == 1
    assert result["ired"]["validation_rows"] == 1
    assert result["ired"]["test_rows"] == 1
    assert "test_structural_scoreability_fraction" in result["ired"]
