from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]

SCALABLE_PATH = ROOT / "v9_trpb_blind" / "scalable_consensus.py"
SPEC = importlib.util.spec_from_file_location("scalable_consensus", SCALABLE_PATH)
SCALABLE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules["scalable_consensus"] = SCALABLE
SPEC.loader.exec_module(SCALABLE)

LEGACY_PATH = ROOT / "v9_pair_transfer_dev" / "run_v9_simple_ablation.py"
LEGACY_SPEC = importlib.util.spec_from_file_location(
    "run_v9_simple_ablation",
    LEGACY_PATH,
)
LEGACY = importlib.util.module_from_spec(LEGACY_SPEC)
assert LEGACY_SPEC.loader is not None
sys.modules["run_v9_simple_ablation"] = LEGACY
LEGACY_SPEC.loader.exec_module(LEGACY)


def _compare(b2, context, ids):
    old_score, old_layer = LEGACY.pareto_consensus_score(
        np.asarray(b2, dtype=float),
        np.asarray(context, dtype=float),
        list(ids),
    )
    new_score, new_layer, _, _ = SCALABLE.pareto_consensus_score_scalable(
        np.asarray(b2, dtype=float),
        np.asarray(context, dtype=float),
        list(ids),
    )
    assert np.array_equal(new_layer, old_layer)
    assert np.array_equal(new_score, old_score)


def test_scalable_matches_legacy_on_random_scores():
    rng = np.random.default_rng(161)
    for n in (5, 17, 61):
        b2 = rng.normal(size=n)
        context = rng.normal(size=n)
        ids = [f"id-{i:03d}" for i in range(n)]
        _compare(b2, context, ids)


def test_scalable_matches_legacy_with_ties():
    b2 = [3, 3, 2, 2, 2, 1, 0, 0]
    context = [3, 2, 3, 2, 2, 4, 0, 0]
    ids = ["z", "a", "b", "c", "d", "e", "f", "g"]
    _compare(b2, context, ids)


def test_scalable_identical_points_share_layer():
    b2 = np.asarray([4.0, 4.0, 3.0])
    context = np.asarray([4.0, 4.0, 3.0])
    score, layer, _, _ = SCALABLE.pareto_consensus_score_scalable(
        b2,
        context,
        ["b", "a", "c"],
    )
    assert layer[0] == layer[1] == 0
    assert layer[2] == 1
    assert score[1] > score[0]
