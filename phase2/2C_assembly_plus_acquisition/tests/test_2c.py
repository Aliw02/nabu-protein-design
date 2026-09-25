from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd


HERE = Path(__file__).resolve().parents[1]
PHASE2 = HERE.parent
PHASE2A0 = PHASE2 / "2A_closed_loop_acquisition" / "2A0_campaign_simulator"
PHASE2A1 = PHASE2 / "2A_closed_loop_acquisition" / "2A1_baselines"
PHASE2B = PHASE2 / "2B_candidate_assembly"
for path in (HERE, PHASE2A0, PHASE2A1, PHASE2B):
    sys.path.insert(0, str(path))

from gb1_adapter import (  # noqa: E402
    GB1_MUTABLE_POSITIONS,
    GB1_REFERENCE_SEQUENCE,
    GB1_WT_GENOTYPE,
    adapt_gb1_frame,
    gb1_reference,
    genotype_to_mutation_set,
)


def test_gb1_reference_matches_wild_type_sites():
    reference = gb1_reference()
    assert len(GB1_REFERENCE_SEQUENCE) == 56
    observed = "".join(
        reference.residue_at(position)
        for position in GB1_MUTABLE_POSITIONS
    )
    assert observed == GB1_WT_GENOTYPE == "VDGV"


def test_gb1_genotype_encoding_uses_true_positions():
    assert genotype_to_mutation_set("ADGV") == ("V39A",)
    assert genotype_to_mutation_set("VAGV") == ("D40A",)
    assert genotype_to_mutation_set("VDAV") == ("G41A",)
    assert genotype_to_mutation_set("VDGA") == ("V54A",)


def test_gb1_wt_is_empty_and_excluded_from_adapter():
    raw = pd.DataFrame(
        {
            "variant": ["VDGV", "ADGV"],
            "fitness": [1.0, 2.0],
        }
    )
    adapted = adapt_gb1_frame(raw)
    assert adapted["candidate_id"].tolist() == ["V39A"]


def test_adapter_is_label_independent_for_identity():
    left = pd.DataFrame(
        {
            "variant": ["ADGV", "VAGV", "VDAV"],
            "fitness": [1.0, 2.0, 3.0],
        }
    )
    right = left.copy()
    right["fitness"] = [-100.0, -200.0, -300.0]

    left_ids = adapt_gb1_frame(left)["candidate_id"].tolist()
    right_ids = adapt_gb1_frame(right)["candidate_id"].tolist()
    assert left_ids == right_ids
