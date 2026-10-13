"""Frozen same-information predictors; no external labels are accepted here."""
from __future__ import annotations

import importlib.util
from itertools import combinations, product
from time import perf_counter

import numpy as np
from scipy.linalg import solve
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR

from common import ALPHABET, HERE, folds_for

LEGACY_PATH = HERE.parent / "v8_1_crossfit_higher_order_dev/run_v8_1_crossfit_dev.py"
spec = importlib.util.spec_from_file_location("frozen_nabu_hierarchy", LEGACY_PATH)
legacy = importlib.util.module_from_spec(spec)
spec.loader.exec_module(legacy)


def one_hot(x):
    return np.eye(len(ALPHABET), dtype=float)[x].reshape(len(x), -1)


class CategoricalANOVA:
    """Joint ridge in centered categorical main and optional pair features."""

    def __init__(self, alpha=1.0, pairs=True):
        self.alpha, self.pairs = alpha, pairs

    def components(self, x, z):
        main = np.zeros((len(x), len(z)))
        squared = np.zeros_like(main)
        for position in range(x.shape[1]):
            p = self.frequencies[position]
            site = ((x[:, position, None] == z[None, :, position]).astype(float)
                    - p[x[:, position, None]] - p[z[None, :, position]]
                    + np.dot(p, p))
            main += site
            squared += site * site
        return main, (main * main - squared) / 2.0

    def kernel(self, x, z):
        main, pair = self.components(x, z)
        kernel = main / self.main_scale if self.main_scale > 1e-12 else np.zeros_like(main)
        if self.pairs and self.pair_scale > 1e-12:
            kernel = kernel + pair / self.pair_scale
        return kernel

    def fit(self, x, y):
        self.x = np.array(x, copy=True)
        self.frequencies = np.stack([
            np.bincount(x[:, j], minlength=len(ALPHABET)) / len(x)
            for j in range(x.shape[1])
        ])
        main, pair = self.components(x, x)
        self.main_scale = float(np.diag(main).mean())
        self.pair_scale = float(np.diag(pair).mean())
        kernel = self.kernel(x, x)
        self.column_mean = kernel.mean(axis=0)
        self.grand_mean = float(kernel.mean())
        centered = kernel - self.column_mean[None, :] - self.column_mean[:, None] + self.grand_mean
        self.y_mean = float(np.mean(y))
        self.coef = solve(centered + self.alpha * np.eye(len(x)), y - self.y_mean,
                          assume_a="pos", check_finite=True)
        return self

    def predict(self, x):
        kernel = self.kernel(x, self.x)
        centered = kernel - kernel.mean(axis=1, keepdims=True) - self.column_mean + self.grand_mean
        return self.y_mean + centered @ self.coef


class LegacyMemory:
    """Only a representation adapter around the unchanged B2/B3 equations."""

    def __init__(self, order):
        self.order = order

    def tokens(self, x):
        return [tuple(f"p{j:03d}:{ALPHABET[self.reference[j]]}>{ALPHABET[state]}"
                      for j, state in enumerate(row) if state != self.reference[j]) for row in x]

    def fit(self, x, y):
        self.reference = np.array([np.bincount(x[:, j], minlength=len(ALPHABET)).argmax()
                                   for j in range(x.shape[1])])
        self.base = legacy.fit_base(self.tokens(x), y)
        return self

    def predict(self, x):
        return np.array([legacy.score_base(tokens, self.base)[self.order - 2] for tokens in self.tokens(x)])

    def support(self, x):
        records = []
        for tokens in self.tokens(x):
            pairs = list(combinations(tokens, 2))
            records.append({"mutation_tokens": len(tokens),
                            "supported_main_tokens": sum(t in self.base["main"] for t in tokens),
                            "possible_pairs": len(pairs),
                            "supported_pairs": sum(t in self.base["pair"] for t in pairs)})
        return records


class NearestNeighbor:
    def fit(self, x, y):
        self.x, self.y = np.array(x, copy=True), np.array(y, copy=True)
        return self

    def predict(self, x):
        distances = (x[:, None, :] != self.x[None, :, :]).sum(axis=2)
        return np.array([self.y[row == row.min()].mean() for row in distances])


class MeanPredictor:
    def fit(self, x, y):
        self.mean = float(np.mean(y))
        return self

    def predict(self, x):
        return np.full(len(x), self.mean)


def features(method, x, lengths):
    if method in ("RBF_SVR", "RANDOM_FOREST"):
        return one_hot(x)
    return lengths if method == "LENGTH_ONLY" else x


def parameter_grid(method, config):
    if method in ("NABU_JOINT_PAIRS", "JOINT_ADDITIVE", "SHUFFLED_JOINT_PAIRS", "LENGTH_ONLY"):
        return [{"alpha": a} for a in config["ridge_alphas"]]
    if method == "RBF_SVR":
        return [{"C": c, "gamma_scale": g} for c, g in product(config["svr_C"], config["svr_gamma"])]
    if method == "RANDOM_FOREST":
        return [{"min_samples_leaf": l, "max_features": f}
                for l, f in product(config["rf_min_samples_leaf"], config["rf_max_features"])]
    return [{}]


def create_model(method, parameters, width, config, seed):
    if method in ("NABU_JOINT_PAIRS", "SHUFFLED_JOINT_PAIRS", "JOINT_ADDITIVE"):
        return CategoricalANOVA(**parameters, pairs=method != "JOINT_ADDITIVE")
    if method in ("LEGACY_B2", "LEGACY_B3"):
        return LegacyMemory(int(method[-1]))
    if method == "LENGTH_ONLY":
        return make_pipeline(StandardScaler(), Ridge(**parameters))
    if method == "RBF_SVR":
        return SVR(C=parameters["C"], gamma=parameters["gamma_scale"] / width, epsilon=0.1)
    if method == "RANDOM_FOREST":
        return RandomForestRegressor(n_estimators=config["rf_trees"], random_state=seed,
                                     n_jobs=1, **parameters)
    if method == "NEAREST_NEIGHBOR":
        return NearestNeighbor()
    if method == "TRAIN_MEAN":
        return MeanPredictor()
    raise ValueError(method)


def fit_predict(method, x, y, lengths, groups, query, query_lengths, config, seed):
    """Nested model selection sees only this training partition's targets."""
    started = perf_counter()
    targets = np.array(y, copy=True)
    if method == "SHUFFLED_JOINT_PAIRS":
        targets = np.random.default_rng(seed).permutation(targets)
    transformed = features(method, x, lengths)
    query_features = features(method, query, query_lengths)
    grid, traces = parameter_grid(method, config), []
    if len(grid) > 1:
        splits = folds_for(groups, config["inner_folds"])
        for parameters in grid:
            predictions = np.full(len(y), np.nan)
            for fit, hold in splits:
                assert not set(groups[fit]).intersection(groups[hold])
                model = create_model(method, parameters, x.shape[1], config, seed)
                model.fit(transformed[fit], targets[fit])
                predictions[hold] = model.predict(transformed[hold])
            assert np.isfinite(predictions).all()
            traces.append({"parameters": parameters,
                           "inner_pooled_mse": float(np.mean((predictions - targets) ** 2))})
        # Grid order is the preregistered deterministic tie breaker.
        selected = min(range(len(traces)), key=lambda i: (traces[i]["inner_pooled_mse"], i))
    else:
        selected = 0
    chosen = grid[selected]
    model = create_model(method, chosen, x.shape[1], config, seed)
    model.fit(transformed, targets)
    predictions = model.predict(query_features)
    assert np.isfinite(predictions).all()
    diagnostics = {"selected_parameters": chosen, "tuning_trace": traces,
                   "fit_and_predict_seconds": perf_counter() - started,
                   "training_records": len(y), "query_records": len(query)}
    support = model.support(query) if isinstance(model, LegacyMemory) else None
    return predictions, diagnostics, support
