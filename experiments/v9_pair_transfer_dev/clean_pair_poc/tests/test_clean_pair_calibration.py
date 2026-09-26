from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import numpy as np
import pandas as pd


RUNNER = (
    Path(__file__).resolve().parents[1]
    / "run_clean_pair_calibration.py"
)
SPEC = importlib.util.spec_from_file_location(
    "run_clean_pair_calibration",
    RUNNER,
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules["run_clean_pair_calibration"] = MODULE
SPEC.loader.exec_module(MODULE)


def test_inner_fold_is_deterministic_and_outer_specific():
    candidate = "ACDEFGHIK"
    first = MODULE.inner_fold(candidate, 0)
    second = MODULE.inner_fold(candidate, 0)
    other_outer = MODULE.inner_fold(candidate, 1)

    assert first == second
    assert 0 <= first < MODULE.N_FOLDS
    assert 0 <= other_outer < MODULE.N_FOLDS


def test_affine_calibration_recovers_exact_linear_map():
    raw = np.asarray([-2.0, -1.0, 0.0, 1.0, 2.0], dtype=float)
    target = 0.25 + 0.4 * raw

    model = MODULE.fit_affine(raw, target)

    assert np.isclose(model["alpha"], 0.25)
    assert np.isclose(model["beta"], 0.4)
    assert model["rank"] == 2


def test_affine_calibration_shrinks_overdispersed_prediction():
    raw = np.asarray([-4.0, -2.0, 0.0, 2.0, 4.0], dtype=float)
    target = np.asarray([-1.0, -0.5, 0.0, 0.5, 1.0], dtype=float)

    model = MODULE.fit_affine(raw, target)

    assert 0.0 < model["beta"] < 1.0
    calibrated = model["alpha"] + model["beta"] * raw
    assert np.sqrt(np.mean((calibrated - target) ** 2)) < np.sqrt(
        np.mean((raw - target) ** 2)
    )


def test_inner_oof_calibration_uses_only_outer_train_rows():
    rows = []
    mutations = ["A1C", "B2D", "C3E", "D4F", "E5G", "F6H"]
    index = 0
    for left_i in range(len(mutations)):
        for right_i in range(left_i + 1, len(mutations)):
            left = mutations[left_i]
            right = mutations[right_i]
            residual = 0.2 + 0.1 * left_i - 0.05 * right_i
            rows.append(
                {
                    "candidate_id": f"candidate_{index}",
                    "left": left,
                    "right": right,
                    "clean_pair_residual": residual,
                    "fold": index % MODULE.N_FOLDS,
                }
            )
            index += 1

    clean = pd.DataFrame(rows)
    outer_fold = 2
    outer_train = clean[clean["fold"] != outer_fold].copy()
    outer_hold_ids = set(
        clean[clean["fold"] == outer_fold]["candidate_id"]
    )

    calibration, inner_oof = MODULE.generate_inner_oof_calibration(
        outer_train,
        outer_fold,
    )

    assert len(inner_oof) > 0
    assert set(inner_oof["candidate_id"]).isdisjoint(outer_hold_ids)
    assert calibration["row_count"] == len(inner_oof)
    assert np.isfinite(calibration["alpha"])
    assert np.isfinite(calibration["beta"])
