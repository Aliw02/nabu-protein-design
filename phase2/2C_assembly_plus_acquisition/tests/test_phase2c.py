from __future__ import annotations

from pathlib import Path
import sys
from types import SimpleNamespace

import pandas as pd


HERE = Path(__file__).resolve().parents[1]
PHASE2 = HERE.parent
PHASE2A0 = (
    PHASE2
    / "2A_closed_loop_acquisition"
    / "2A0_campaign_simulator"
)
PHASE2A1 = (
    PHASE2
    / "2A_closed_loop_acquisition"
    / "2A1_baselines"
)
PHASE2B = PHASE2 / "2B_candidate_assembly"

for path in (
    HERE,
    PHASE2A0,
    PHASE2A1,
    PHASE2B,
):
    sys.path.insert(0, str(path))

from four_site_codec import (  # noqa: E402
    FourSiteCodec,
    deterministic_identity_subpool,
)
from maturity_gate import proposal_is_mature  # noqa: E402
from run_multilandscape_benchmark import (  # noqa: E402
    LANDSCAPE_SPECS,
    _joint_degradation,
    _landscape_win,
)


def test_four_site_codec_roundtrip():
    codec = FourSiteCodec("VDGV")
    mutation_set = codec.encode("ADAV")
    assert mutation_set == (
        "V1A",
        "G3A",
    )
    assert codec.decode(mutation_set) == "ADAV"


def test_reference_genotype_encodes_to_empty_set():
    codec = FourSiteCodec("VFVS")
    assert codec.encode("VFVS") == ()
    assert codec.decode(()) == "VFVS"


def test_identity_subpool_is_truth_independent():
    codec = FourSiteCodec("VDGV")
    genotypes = [
        "VDGV",
        "ADGV",
        "CDGV",
        "DDGV",
        "EDGV",
        "FDGV",
        "GDGV",
        "HDGV",
        "IDGV",
        "KDGV",
    ]
    frame = pd.DataFrame(
        {
            "candidate_id": genotypes,
            "mutation_set": [
                codec.encode(value)
                for value in genotypes
            ],
            "crossfit_fold": [
                index % 5
                for index in range(len(genotypes))
            ],
            "DMS_score": [
                float(index)
                for index in range(len(genotypes))
            ],
        }
    )
    altered = frame.copy()
    altered["DMS_score"] = list(
        reversed(
            altered["DMS_score"].tolist()
        )
    )

    left = deterministic_identity_subpool(
        frame,
        pool_size=6,
    )
    right = deterministic_identity_subpool(
        altered,
        pool_size=6,
    )
    assert left["candidate_id"].tolist() == right[
        "candidate_id"
    ].tolist()


def test_maturity_requires_repeated_main_and_pair_support():
    model = SimpleNamespace(
        model={
            "base": {
                "main": {
                    "V1A": {"support": 2},
                    "D2A": {"support": 2},
                    "G3A": {"support": 2},
                },
                "pair": {
                    ("V1A", "D2A"): {
                        "support": 2
                    },
                    ("V1A", "G3A"): {
                        "support": 2
                    },
                    ("D2A", "G3A"): {
                        "support": 2
                    },
                },
            }
        }
    )
    assert proposal_is_mature(
        ("V1A", "D2A", "G3A"),
        model,
        min_support=2,
    )

    model.model["base"]["pair"][
        ("D2A", "G3A")
    ]["support"] = 1
    assert not proposal_is_mature(
        ("V1A", "D2A", "G3A"),
        model,
        min_support=2,
    )


def _summary(auc, regret, hits, assembly_rounds):
    return {
        "post_bootstrap_normalized_discovery_auc": auc,
        "assembly_active_rounds": assembly_rounds,
        "final": {
            "best_so_far_regret": regret,
            "cumulative_top1_hits": hits,
        },
    }


def test_landscape_win_rule_is_preregistered_and_strict():
    acquisition = _summary(
        0.70,
        0.20,
        3,
        0,
    )
    combined = _summary(
        0.75,
        0.20,
        3,
        2,
    )
    assert _landscape_win(
        acquisition,
        combined,
    )

    combined["post_bootstrap_normalized_discovery_auc"] = 0.70
    assert not _landscape_win(
        acquisition,
        combined,
    )


def test_joint_degradation_requires_auc_and_regret_worse():
    acquisition = _summary(
        0.70,
        0.20,
        3,
        0,
    )
    combined = _summary(
        0.69,
        0.30,
        4,
        2,
    )
    assert _joint_degradation(
        acquisition,
        combined,
    )


def test_phase2c_specs_do_not_include_reserved_blind_sets():
    assert set(LANDSCAPE_SPECS) == {
        "GB1",
        "TRPB",
        "PHOQ",
    }
    text = repr(LANDSCAPE_SPECS).upper()
    assert "AAV" not in text
    assert "IRED" not in text
