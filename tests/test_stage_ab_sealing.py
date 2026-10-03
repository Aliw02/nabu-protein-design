from __future__ import annotations

import hashlib
from itertools import combinations

import numpy as np
import pandas as pd

from nabu_protein.higher_order import fnv1a32
from nabu_protein.metrics import evaluate_v83
from nabu_protein.v83 import NabuV83Model


def build_synthetic_landscape():
    mutations = (
        "A1C",
        "B2D",
        "C3E",
        "D4F",
        "E5G",
        "F6H",
        "G7I",
        "H8K",
    )
    mutation_sets = (
        list(combinations(mutations, 3))
        + list(combinations(mutations, 4))
    )
    rows = []
    for i, mutation_set in enumerate(mutation_sets):
        candidate_id = ":".join(mutation_set)
        truth = (
            0.15 * len(mutation_set)
            + 0.07 * sum(int(token[1:-1]) for token in mutation_set)
            + 0.013 * ((i * 11) % 17)
        )
        rows.append(
            {
                "candidate_id": candidate_id,
                "mutation_set": mutation_set,
                "truth": truth,
                "bucket": fnv1a32(candidate_id) % 10,
            }
        )
    return pd.DataFrame(rows)


def freeze_manifest(frame):
    ranked = frame.sort_values(
        ["V8_3_ADAPTIVE_ROUTER", "candidate_id"],
        ascending=[False, True],
    )
    payload = "\n".join(
        f"{row.candidate_id},{row.V8_3_ADAPTIVE_ROUTER:.17g}"
        for row in ranked.itertuples()
    ).encode("utf-8")
    return {
        "candidate_ids": ranked["candidate_id"].tolist(),
        "sha256": hashlib.sha256(payload).hexdigest(),
    }


def test_stage_a_predictions_are_frozen_before_stage_b_truth():
    landscape = build_synthetic_landscape()
    visible = landscape[landscape["bucket"] <= 6].copy()
    hidden = landscape[landscape["bucket"] >= 7].copy()

    model = NabuV83Model().fit(
        visible["mutation_set"].tolist(),
        visible["truth"].to_numpy(dtype=float),
        visible["candidate_id"].tolist(),
    )

    stage_a = model.score_candidates(
        hidden["mutation_set"].tolist(),
        hidden["candidate_id"].tolist(),
    )
    stage_a = stage_a[stage_a["scoreable"]].copy()
    frozen_before_reveal = freeze_manifest(stage_a)

    # Stage B reveal 1: the actual synthetic truth.
    reveal_1 = stage_a.merge(
        hidden[["candidate_id", "truth"]],
        on="candidate_id",
        validate="one_to_one",
    )
    metrics_1 = evaluate_v83(
        reveal_1,
        "V8_3_ADAPTIVE_ROUTER",
        "truth",
    )

    # Stage B reveal 2: deliberately different hidden truth. Predictions are
    # not recomputed; only evaluation changes.
    reveal_2 = reveal_1.copy()
    reveal_2["truth"] = -reveal_2["truth"].to_numpy(dtype=float)
    metrics_2 = evaluate_v83(
        reveal_2,
        "V8_3_ADAPTIVE_ROUTER",
        "truth",
    )

    frozen_after_reveal = freeze_manifest(stage_a)

    assert frozen_after_reveal == frozen_before_reveal
    assert metrics_1 != metrics_2
