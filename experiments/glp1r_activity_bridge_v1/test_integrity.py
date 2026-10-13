"""Pre-outcome tests for mathematical equivalence and split isolation."""
from __future__ import annotations

import os
for variable in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[variable] = "1"

from datetime import datetime, timezone
from itertools import combinations

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge

from common import ALPHABET, HERE, codes_from, folds_for, read_json, sequence_groups, write_json
from metrics import metrics, paired_bootstrap
from models import CategoricalANOVA, LegacyMemory, fit_predict, legacy


def test_kernel_matches_explicit_features_and_primal_ridge():
    rng = np.random.default_rng(714)
    x = rng.integers(0, 3, size=(25, 4))
    query = rng.integers(0, 4, size=(11, 4))
    y = rng.normal(size=len(x))
    for use_pairs in [False, True]:
        model = CategoricalANOVA(alpha=0.3, pairs=use_pairs).fit(x, y)
        def explicit(values):
            site = [np.eye(len(ALPHABET))[values[:, j]] - model.frequencies[j] for j in range(x.shape[1])]
            blocks = [np.concatenate(site, axis=1) / np.sqrt(model.main_scale)]
            if use_pairs:
                blocks.append(np.concatenate([(site[i][:, :, None] * site[j][:, None, :]).reshape(len(values), -1)
                                               for i, j in combinations(range(x.shape[1]), 2)], axis=1) / np.sqrt(model.pair_scale))
            return np.concatenate(blocks, axis=1)
        train_features, test_features = explicit(x), explicit(query)
        np.testing.assert_allclose(model.kernel(x, x), train_features @ train_features.T, atol=1e-12)
        np.testing.assert_allclose(model.kernel(query, x), test_features @ train_features.T, atol=1e-12)
        assert np.linalg.eigvalsh(model.kernel(x, x)).min() > -1e-10
        primal = Ridge(alpha=0.3, fit_intercept=True, solver="svd").fit(train_features, y)
        np.testing.assert_allclose(model.predict(query), primal.predict(test_features), atol=1e-10)
        np.testing.assert_allclose(model.predict(query), np.array([model.predict(row[None, :])[0] for row in query]), atol=1e-12)
        shifted = CategoricalANOVA(alpha=0.3, pairs=use_pairs).fit(x, y + 7)
        np.testing.assert_allclose(shifted.predict(query), model.predict(query) + 7, atol=1e-11)
    constant = np.zeros((9, 4), dtype=int)
    result = CategoricalANOVA().fit(constant, np.arange(9.)).predict(query)
    np.testing.assert_allclose(result, 4)


def test_legacy_uses_canonical_equations_without_clipping():
    x = np.array([[0, 0, 0], [0, 1, 0], [1, 0, 0], [1, 1, 0], [1, 1, 1]])
    y = np.array([0., 1., 3., 8., -2.])
    query = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2]])
    for order in [2, 3]:
        model = LegacyMemory(order).fit(x, y)
        reference = legacy.fit_base(model.tokens(x), y)
        expected = np.array([legacy.score_base(tokens, reference)[order - 2] for tokens in model.tokens(query)])
        np.testing.assert_array_equal(model.predict(query), expected)
        assert model.support(query)[-1]["supported_main_tokens"] == 0
        assert model.support(query)[-1]["supported_pairs"] == 0


def test_real_identity_group_isolation_without_external_targets():
    dev = pd.read_csv(HERE / "data/prepared/development.csv")
    ext = pd.read_csv(HERE / "data/prepared/external_identities.csv")
    config = read_json(HERE / "protocol.json")
    groups, distances = sequence_groups(dev.sequence.tolist(), config["cluster_edit_distance"])
    np.testing.assert_array_equal(groups, dev.sequence_group)
    seen = np.zeros(len(dev), dtype=int)
    for fold, (fit, hold) in enumerate(folds_for(groups, config["outer_folds"])):
        assert (dev.outer_fold.iloc[hold] == fold).all()
        assert distances[np.ix_(fit, hold)].min() > 2
        seen[hold] += 1
        for inner_fit, inner_hold in folds_for(groups[fit], config["inner_folds"]):
            assert distances[np.ix_(fit[inner_fit], fit[inner_hold])].min() > 2
    assert (seen == 1).all()
    for frame in [dev, ext]:
        assert frame.encoded_sequence.str.replace("-", "", regex=False).eq(frame.sequence).all()
        assert np.isfinite(codes_from(frame)).all()
    assert not set(dev.sequence).intersection(ext.sequence)
    assert "target_pEC50" not in ext


def test_selection_does_not_depend_on_query_composition():
    rng = np.random.default_rng(217)
    x = rng.integers(0, 4, size=(24, 4))
    y = rng.normal(size=24)
    groups, lengths = np.repeat(np.arange(8), 3), np.ones((24, 3))
    query, query_lengths = rng.integers(0, 5, size=(5, 4)), np.ones((5, 3))
    config = read_json(HERE / "protocol.json")
    for method in ["NABU_JOINT_PAIRS", "JOINT_ADDITIVE", "SHUFFLED_JOINT_PAIRS"]:
        first, trace, _ = fit_predict(method, x, y, lengths, groups, query, query_lengths, config, 17)
        second, second_trace, _ = fit_predict(method, x, y, lengths, groups, query[:1], query_lengths[:1], config, 17)
        np.testing.assert_allclose(first[:1], second, atol=1e-12)
        assert trace["tuning_trace"] == second_trace["tuning_trace"]
        assert trace["selected_parameters"] == second_trace["selected_parameters"]


def test_metrics_ties_and_paired_bootstrap_direction():
    y = np.arange(10., dtype=float)
    ideal = metrics(y, y, [str(i) for i in range(10)])
    assert ideal["rmse"] == 0 and ideal["top5_recall"] == 1
    assert abs(ideal["spearman"] - 1) < 1e-12
    tied = metrics(y, np.ones(10), [f"id{i:02d}" for i in range(10)])
    assert tied["spearman"] is None and tied["top5_ids"] == [f"id{i:02d}" for i in range(5)]
    values = paired_bootstrap(y, y, y + 2, np.repeat(np.arange(5), 2), 100, 3)
    assert values["ci95_low"] == values["ci95_high"] == -4 and values["upper_bound_negative"]


if __name__ == "__main__":
    tests = [value for name, value in sorted(globals().items()) if name.startswith("test_")]
    for test in tests:
        test()
        print(f"PASS {test.__name__}")
    report = {"passed": True, "tests": [test.__name__ for test in tests],
              "tested_utc": datetime.now(timezone.utc).isoformat(),
              "synthetic_model_fits_only": True, "external_outcomes_loaded": False}
    # Keep the initial test record immutable after the first complete freeze.
    if not (HERE / "results/PREDICTION_FREEZE.json").exists():
        write_json(HERE / "INTEGRITY_TESTS.json", report)
