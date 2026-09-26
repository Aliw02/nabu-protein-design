from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import numpy as np


RUNNER = Path(__file__).resolve().parents[1] / "run_v9_simple_ablation.py"
SPEC = importlib.util.spec_from_file_location("run_v9_simple_ablation", RUNNER)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules["run_v9_simple_ablation"] = MODULE
SPEC.loader.exec_module(MODULE)


def test_pareto_consensus_prefers_jointly_dominant_candidate():
    b2 = np.asarray([3.0, 2.0, 1.0], dtype=float)
    context = np.asarray([3.0, 1.0, 2.0], dtype=float)
    ids = ["a", "b", "c"]

    score, layer = MODULE.pareto_consensus_score(b2, context, ids)

    assert layer[0] == 0
    assert score[0] == np.max(score)


def test_pareto_consensus_is_deterministic():
    b2 = np.asarray([1.0, 2.0, 3.0, 2.0], dtype=float)
    context = np.asarray([3.0, 2.0, 1.0, 2.0], dtype=float)
    ids = ["a", "b", "c", "d"]

    first_score, first_layer = MODULE.pareto_consensus_score(
        b2, context, ids
    )
    second_score, second_layer = MODULE.pareto_consensus_score(
        b2, context, ids
    )

    assert np.array_equal(first_score, second_score)
    assert np.array_equal(first_layer, second_layer)


def test_contact_neighbor_channel_uses_contact_weighted_context():
    batch = 1
    length = 3
    hidden = 320

    reference_hidden = np.zeros((length, hidden), dtype=np.float32)
    residue_hidden = np.zeros((batch, length, hidden), dtype=np.float32)
    residue_hidden[0, 2, :] = 4.0

    mutation_mask = np.asarray([[0, 1, 0]], dtype=np.float32)
    contacts = np.zeros((batch, length, length), dtype=np.float32)
    contacts[0, 1, 2] = 1.0

    feature = MODULE.contact_neighbor_features(
        residue_hidden,
        reference_hidden,
        mutation_mask,
        contacts,
    )

    assert feature.shape == (1, MODULE.CONTACT_FEATURE_DIM)
    assert np.allclose(feature[0], 4.0)


def test_contact_neighbor_wt_is_zero():
    length = 2
    hidden = 320

    reference_hidden = np.zeros((length, hidden), dtype=np.float32)
    residue_hidden = np.ones((1, length, hidden), dtype=np.float32)
    mutation_mask = np.zeros((1, length), dtype=np.float32)
    contacts = np.ones((1, length, length), dtype=np.float32)

    feature = MODULE.contact_neighbor_features(
        residue_hidden,
        reference_hidden,
        mutation_mask,
        contacts,
    )

    assert np.allclose(feature, 0.0)
