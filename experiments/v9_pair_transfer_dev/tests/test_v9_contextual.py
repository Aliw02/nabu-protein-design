from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import numpy as np
import pandas as pd


RUNNER = Path(__file__).resolve().parents[1] / "run_v9_contextual.py"
SPEC = importlib.util.spec_from_file_location("run_v9_contextual", RUNNER)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules["run_v9_contextual"] = MODULE
SPEC.loader.exec_module(MODULE)


def test_contextual_feature_shape_and_wt_zero_mutation_channel():
    batch = 2
    length = 4
    hidden = 320

    reference_hidden = np.zeros((length, hidden), dtype=np.float32)
    residue_hidden = np.zeros((batch, length, hidden), dtype=np.float32)
    residue_hidden[1, 2, :] = 2.0

    mutation_mask = np.asarray(
        [
            [0, 0, 0, 0],
            [0, 0, 1, 0],
        ],
        dtype=np.float32,
    )

    feature = MODULE.contextual_features_from_hidden(
        residue_hidden,
        reference_hidden,
        mutation_mask,
    )

    assert feature.shape == (batch, MODULE.FEATURE_DIM)
    assert np.allclose(feature[0, 320:], 0.0)
    assert np.allclose(feature[1, 320:], 2.0)


def test_same_mutation_identity_can_have_different_context_feature():
    length = 3
    hidden = 320
    reference_hidden = np.zeros((length, hidden), dtype=np.float32)

    residue_hidden = np.zeros((2, length, hidden), dtype=np.float32)
    residue_hidden[0, 1, :] = 1.0
    residue_hidden[1, 1, :] = 3.0

    mutation_mask = np.asarray(
        [
            [0, 1, 0],
            [0, 1, 0],
        ],
        dtype=np.float32,
    )

    feature = MODULE.contextual_features_from_hidden(
        residue_hidden,
        reference_hidden,
        mutation_mask,
    )

    assert not np.allclose(feature[0], feature[1])
    assert np.allclose(feature[0, 320:], 1.0)
    assert np.allclose(feature[1, 320:], 3.0)


def test_clean_additive_baseline_uses_b2_only_when_singleton_missing():
    mutations = [
        ("A1C", "B2D"),
        ("A1C", "X9Y"),
    ]
    baseline, supported = MODULE.build_clean_additive_baseline(
        mutations,
        wt_target=1.0,
        clean_main={
            "A1C": 0.4,
            "B2D": -0.2,
        },
        b2_prediction=np.asarray([9.0, 3.25], dtype=float),
    )

    assert supported.tolist() == [True, False]
    assert np.isclose(baseline[0], 1.2)
    assert np.isclose(baseline[1], 3.25)


def test_residual_training_table_anchors_wt_and_single_to_zero():
    fit = pd.DataFrame(
        {
            "sequence": ["AAA", "CAA", "ACA", "CCA"],
            "target": [1.0, 1.5, 0.8, 1.7],
        }
    )
    clean = pd.DataFrame(
        {
            "candidate_id": ["CCA"],
            "clean_pair_residual": [0.4],
        }
    )

    table = MODULE.build_residual_training_table(
        fit,
        reference="AAA",
        clean_pairs=clean,
    )

    wt = table[table["evidence_class"] == "WT"]
    singles = table[table["evidence_class"] == "SINGLE"]
    doubles = table[table["evidence_class"] == "CLEAN_DOUBLE"]

    assert len(wt) == 1
    assert len(singles) == 2
    assert len(doubles) == 1
    assert np.allclose(wt["interaction_residual"], 0.0)
    assert np.allclose(singles["interaction_residual"], 0.0)
    assert np.isclose(doubles.iloc[0]["interaction_residual"], 0.4)
