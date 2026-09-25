from __future__ import annotations

from itertools import combinations
from pathlib import Path
import sys

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from campaign import (  # noqa: E402
    CampaignSimulator,
    IdentityHashPolicy,
    VirtualAssayOracle,
    bootstrap_seed_ids,
    identity_pool_from_frame,
    selection_sha256,
    truth_frame_from_frame,
)


def synthetic_landscape() -> pd.DataFrame:
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
        list(combinations(mutations, 2))
        + list(combinations(mutations, 3))
    )

    rows = []
    for index, mutation_set in enumerate(mutation_sets):
        candidate_id = ":".join(mutation_set)
        label = (
            0.31 * len(mutation_set)
            + 0.071 * sum(int(token[1:-1]) for token in mutation_set)
            + 0.019 * ((index * 13) % 17)
        )
        rows.append(
            {
                "candidate_id": candidate_id,
                "mutant": candidate_id,
                "DMS_score": float(label),
            }
        )
    return pd.DataFrame(rows)


class TopVisibleModelPolicy:
    """Uses only candidate scores present in PolicyView, never oracle truth."""

    def select(self, view, batch_size):
        work = view.candidates.copy()
        work["_rank_score"] = work["V8_3_ADAPTIVE_ROUTER"].fillna(-np.inf)
        ranked = work.sort_values(
            ["_rank_score", "candidate_id"],
            ascending=[False, True],
        )
        return ranked.head(batch_size)["candidate_id"].astype(str).tolist()


def make_simulator(frame: pd.DataFrame, out=None) -> CampaignSimulator:
    identities = identity_pool_from_frame(frame)
    truth = truth_frame_from_frame(frame)
    return CampaignSimulator(
        identities,
        VirtualAssayOracle(truth),
        output_dir=out,
    )


def test_identity_only_bootstrap_populates_every_crossfit_fold():
    frame = synthetic_landscape()
    identities = identity_pool_from_frame(frame)

    selected = bootstrap_seed_ids(
        identities,
        target_count=24,
        min_per_fold=1,
    )

    folds = identities.set_index("candidate_id").loc[
        selected, "crossfit_fold"
    ]
    assert len(selected) == 24
    assert set(folds.tolist()) == {0, 1, 2, 3, 4}


def test_policy_view_contains_no_truth_columns():
    simulator = make_simulator(synthetic_landscape())
    simulator.bootstrap(seed_count=24)
    view = simulator.policy_view()

    forbidden = {"assay_value", "DMS_score", "truth", "fitness"}
    assert forbidden.isdisjoint(view.candidates.columns)


def test_oracle_reveals_only_requested_ids():
    frame = synthetic_landscape()
    truth = truth_frame_from_frame(frame)
    oracle = VirtualAssayOracle(truth)

    requested = frame["candidate_id"].head(3).tolist()
    revealed = oracle.reveal(requested)

    assert revealed["candidate_id"].tolist() == requested
    assert len(revealed) == 3


def test_selection_is_frozen_before_reveal(tmp_path):
    simulator = make_simulator(synthetic_landscape(), out=tmp_path)
    simulator.bootstrap(seed_count=24)
    before = simulator.measurements_spent

    frozen = simulator.prepare_selection(
        IdentityHashPolicy(),
        batch_size=5,
    )

    assert simulator.measurements_spent == before
    assert frozen.selection_sha256 == selection_sha256(
        list(frozen.candidate_ids)
    )

    selection_file = (
        tmp_path
        / "rounds"
        / "round_000_selection_before_reveal.json"
    )
    text = selection_file.read_text(encoding="utf-8")
    assert "assay_value" not in text
    assert "DMS_score" not in text

    simulator.reveal_and_update(frozen)
    assert simulator.measurements_spent == before + 5


def test_unselected_truth_changes_cannot_change_pre_reveal_selection():
    original = synthetic_landscape()
    identities = identity_pool_from_frame(original)

    seed_ids = bootstrap_seed_ids(
        identities,
        target_count=24,
        min_per_fold=1,
    )
    seed_set = set(seed_ids)

    altered = original.copy()
    mask = ~altered["candidate_id"].isin(seed_set)
    altered.loc[mask, "DMS_score"] = (
        -1000.0
        - np.arange(int(mask.sum()), dtype=float)
    )

    simulator_a = make_simulator(original)
    simulator_b = make_simulator(altered)

    bootstrap_a = simulator_a.bootstrap(seed_count=24)
    bootstrap_b = simulator_b.bootstrap(seed_count=24)

    assert bootstrap_a[
        "selection_sha256_before_label_reveal"
    ] == bootstrap_b[
        "selection_sha256_before_label_reveal"
    ]
    assert bootstrap_a["revealed_labels_sha256"] == bootstrap_b[
        "revealed_labels_sha256"
    ]

    policy = TopVisibleModelPolicy()
    frozen_a = simulator_a.prepare_selection(policy, batch_size=7)
    frozen_b = simulator_b.prepare_selection(policy, batch_size=7)

    assert frozen_a.candidate_ids == frozen_b.candidate_ids
    assert frozen_a.selection_sha256 == frozen_b.selection_sha256


def test_deterministic_replay_reproduces_round_sequence():
    frame = synthetic_landscape()
    policy = IdentityHashPolicy()

    simulator_a = make_simulator(frame)
    simulator_b = make_simulator(frame)
    simulator_a.bootstrap(seed_count=24)
    simulator_b.bootstrap(seed_count=24)

    for _ in range(3):
        simulator_a.run_round(policy, batch_size=6)
        simulator_b.run_round(policy, batch_size=6)

    ids_a = [
        manifest["candidate_ids"]
        for manifest in simulator_a.round_manifests
    ]
    ids_b = [
        manifest["candidate_ids"]
        for manifest in simulator_b.round_manifests
    ]
    hashes_a = [
        manifest["selection_sha256_before_label_reveal"]
        for manifest in simulator_a.round_manifests
    ]
    hashes_b = [
        manifest["selection_sha256_before_label_reveal"]
        for manifest in simulator_b.round_manifests
    ]

    assert ids_a == ids_b
    assert hashes_a == hashes_b
    assert simulator_a.measurements_spent == simulator_b.measurements_spent
