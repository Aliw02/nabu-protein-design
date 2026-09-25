from __future__ import annotations

import importlib.util
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from nabu_protein.higher_order import (
    N_FOLDS,
    crossfit_fold,
    fit_crossfitted_hierarchy,
    score_hierarchy,
)
from nabu_protein.router import apply_router, oof_elite_diagnostic


ROOT = Path(__file__).resolve().parents[1]


def load_module(name: str, relative_path: str):
    path = ROOT / relative_path
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


REF_V81 = load_module(
    "nabu_ref_v81",
    "experiments/v8_1_crossfit_higher_order_dev/run_v8_1_crossfit_dev.py",
)
REF_V83 = load_module(
    "nabu_ref_v83",
    "experiments/v8_3_multilandscape_validation/run_v8_3_validation.py",
)


def synthetic_visible():
    mutations = (
        "A1C",
        "B2D",
        "C3E",
        "D4F",
        "E5G",
        "F6H",
    )
    base_sets = list(combinations(mutations, 4))[:10]
    base_sets += list(combinations(mutations, 3))[:10]

    candidate_ids = []
    mutation_sets = []
    fold_counts = {fold: 0 for fold in range(N_FOLDS)}

    i = 0
    while min(fold_counts.values()) < 8:
        candidate_id = f"SYN_{i:04d}"
        fold = crossfit_fold(candidate_id)
        if fold_counts[fold] < 8:
            candidate_ids.append(candidate_id)
            mutation_sets.append(base_sets[len(candidate_ids) % len(base_sets)])
            fold_counts[fold] += 1
        i += 1

    labels = np.array(
        [
            0.7 * len(ms)
            + 0.11 * (i % 9)
            + 0.03 * ((i * 7) % 5)
            for i, ms in enumerate(mutation_sets)
        ],
        dtype=float,
    )
    return candidate_ids, mutation_sets, labels


def assert_memory_equal(left, right):
    assert left.keys() == right.keys()
    for key in left:
        assert left[key].keys() == right[key].keys()
        for field in left[key]:
            assert left[key][field] == right[key][field]


def test_b2_b3_b4_b5_matches_frozen_reference_exactly():
    candidate_ids, mutation_sets, labels = synthetic_visible()

    packaged = fit_crossfitted_hierarchy(
        mutation_sets,
        labels,
        candidate_ids,
    )
    frozen = REF_V81.fit_crossfitted_hierarchy(
        mutation_sets,
        labels,
        candidate_ids,
    )

    np.testing.assert_array_equal(packaged["folds"], frozen["folds"])
    for key in (
        "oof_b2",
        "oof_b3",
        "oof_b4",
        "residual3_oof",
        "residual4_oof",
    ):
        np.testing.assert_allclose(
            packaged[key],
            frozen[key],
            rtol=0.0,
            atol=0.0,
        )

    assert packaged["base"]["global_mean"] == frozen["base"]["global_mean"]
    assert_memory_equal(packaged["base"]["main"], frozen["base"]["main"])
    assert_memory_equal(packaged["base"]["pair"], frozen["base"]["pair"])
    assert_memory_equal(packaged["triplet"], frozen["triplet"])
    assert_memory_equal(packaged["quartet"], frozen["quartet"])

    for mutation_set in mutation_sets[:12]:
        assert score_hierarchy(mutation_set, packaged) == (
            REF_V81.score_hierarchy(mutation_set, frozen)
        )


def test_router_matches_frozen_reference_semantics():
    candidate_ids, mutation_sets, labels = synthetic_visible()
    model = fit_crossfitted_hierarchy(
        mutation_sets,
        labels,
        candidate_ids,
    )

    visible = pd.DataFrame({"candidate_id": candidate_ids})
    frozen_elite = REF_V83.oof_elite_diagnostic(
        visible,
        labels,
        model,
        fraction=0.20,
    )
    packaged_elite = oof_elite_diagnostic(
        candidate_ids,
        labels,
        model,
        fraction=0.20,
    )

    assert packaged_elite["region_count"] == frozen_elite["region_count"]
    assert packaged_elite["allow_elite_rerank"] == frozen_elite["allow_elite_rerank"]
    assert packaged_elite["B3"] == frozen_elite["B3"]
    assert (
        packaged_elite["delta_top50_mean_percentile"]
        == frozen_elite["delta_top50_mean_percentile"]
    )
    assert (
        packaged_elite["delta_top50_top1pct_hits"]
        == frozen_elite["delta_top50_top1pct_hits"]
    )

    rows = [score_hierarchy(ms, model) for ms in mutation_sets]
    frame = pd.DataFrame(rows)
    frame["candidate_id"] = candidate_ids

    b3_oof = float(spearmanr(model["oof_b3"], labels).statistic)
    b4_oof = float(spearmanr(model["oof_b4"], labels).statistic)
    if b4_oof - b3_oof > 0:
        expected_mode = "GLOBAL_HIGHER_ORDER"
        expected = frame["B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER"].astype(float)
    elif frozen_elite["allow_elite_rerank"]:
        expected_mode = "RANK_PRESERVING_B3_TOP20_B5_RERANK"
        expected = REF_V83.rank_preserving_elite_rerank(
            frame.copy(),
            fraction=0.20,
        )
    else:
        expected_mode = "B3_PROTECTED_NO_HIGHER_ORDER"
        expected = frame["B3_RAW_PAIR"].astype(float)

    actual_frame = frame.copy()
    decision = apply_router(
        actual_frame,
        model,
        candidate_ids,
        labels,
        "B3_RAW_PAIR",
        "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER",
        "V8_3_ADAPTIVE_ROUTER",
    )

    assert decision["mode"] == expected_mode
    np.testing.assert_allclose(
        actual_frame["V8_3_ADAPTIVE_ROUTER"].to_numpy(dtype=float),
        expected.to_numpy(dtype=float),
        rtol=0.0,
        atol=0.0,
    )
