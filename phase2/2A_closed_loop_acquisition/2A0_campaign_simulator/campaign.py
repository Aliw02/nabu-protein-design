from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
from typing import Protocol

import numpy as np
import pandas as pd

from nabu_protein.higher_order import (
    FROZEN_CORE_COMMIT,
    N_FOLDS,
    crossfit_fold,
    fnv1a32,
    model_diagnostics,
)
from nabu_protein.v83 import NabuV83Model, canonicalize_mutation_set


LABEL_COLUMN = "assay_value"


def parse_mutant_string(value: str) -> tuple[str, ...]:
    text = str(value).strip()
    if not text:
        raise ValueError("Mutation string must not be blank.")
    parts = tuple(part.strip() for part in text.split(":"))
    if any(not part for part in parts):
        raise ValueError(f"Malformed mutation string: {value!r}")
    return canonicalize_mutation_set(parts)


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def selection_sha256(candidate_ids: list[str] | tuple[str, ...]) -> str:
    payload = "\n".join(str(candidate_id) for candidate_id in candidate_ids)
    return sha256_text(payload)


def revealed_labels_sha256(frame: pd.DataFrame) -> str:
    ordered = frame[["candidate_id", LABEL_COLUMN]].copy()
    payload = "\n".join(
        f"{row.candidate_id},{float(row.assay_value):.17g}"
        for row in ordered.itertuples(index=False)
    )
    return sha256_text(payload)


def identity_pool_from_frame(
    frame: pd.DataFrame,
    candidate_column: str = "candidate_id",
    mutation_column: str = "mutant",
) -> pd.DataFrame:
    required = {candidate_column, mutation_column}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Missing identity columns: {sorted(missing)}")

    pool = frame[[candidate_column, mutation_column]].copy()
    pool.columns = ["candidate_id", "mutant"]
    pool["candidate_id"] = pool["candidate_id"].map(str)

    if pool["candidate_id"].map(lambda value: not value.strip()).any():
        raise ValueError("Candidate IDs must not be blank.")
    if pool["candidate_id"].duplicated().any():
        raise ValueError("Candidate IDs must be unique.")

    pool["mutation_set"] = pool["mutant"].map(parse_mutant_string)
    if pool["mutation_set"].duplicated().any():
        raise ValueError(
            "Canonical mutation sets must be unique in the campaign pool."
        )
    pool["crossfit_fold"] = pool["candidate_id"].map(crossfit_fold).astype(int)
    return pool.reset_index(drop=True)


def truth_frame_from_frame(
    frame: pd.DataFrame,
    candidate_column: str = "candidate_id",
    label_column: str = "DMS_score",
) -> pd.DataFrame:
    required = {candidate_column, label_column}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Missing oracle columns: {sorted(missing)}")

    truth = frame[[candidate_column, label_column]].copy()
    truth.columns = ["candidate_id", LABEL_COLUMN]
    truth["candidate_id"] = truth["candidate_id"].map(str)
    truth[LABEL_COLUMN] = truth[LABEL_COLUMN].astype(float)

    if truth["candidate_id"].duplicated().any():
        raise ValueError("Oracle candidate IDs must be unique.")
    if not np.isfinite(truth[LABEL_COLUMN].to_numpy(dtype=float)).all():
        raise ValueError("Oracle labels must contain only finite values.")
    return truth.reset_index(drop=True)


def bootstrap_seed_ids(
    identity_pool: pd.DataFrame,
    target_count: int,
    min_per_fold: int = 1,
    namespace: str = "NABU_PHASE2A0_BOOTSTRAP",
) -> list[str]:
    if target_count < 1:
        raise ValueError("target_count must be positive.")
    if min_per_fold < 1:
        raise ValueError("min_per_fold must be positive.")

    required_count = N_FOLDS * min_per_fold
    target_count = max(int(target_count), required_count)
    if target_count > len(identity_pool):
        raise ValueError(
            f"Bootstrap requires {target_count} candidates but pool has "
            f"{len(identity_pool)}."
        )

    work = identity_pool[["candidate_id", "crossfit_fold"]].copy()
    work["_bootstrap_hash"] = work["candidate_id"].map(
        lambda candidate_id: fnv1a32(
            f"{namespace}|{candidate_id}"
        )
    )
    work = work.sort_values(
        ["_bootstrap_hash", "candidate_id"],
        ascending=[True, True],
    ).reset_index(drop=True)

    selected: list[str] = []
    selected_set: set[str] = set()

    for fold in range(N_FOLDS):
        candidates = work[work["crossfit_fold"] == fold]
        if len(candidates) < min_per_fold:
            raise ValueError(
                f"Fold {fold} has only {len(candidates)} candidates; "
                f"requires {min_per_fold}."
            )
        for candidate_id in candidates.head(min_per_fold)["candidate_id"]:
            candidate_id = str(candidate_id)
            selected.append(candidate_id)
            selected_set.add(candidate_id)

    for candidate_id in work["candidate_id"]:
        candidate_id = str(candidate_id)
        if candidate_id in selected_set:
            continue
        selected.append(candidate_id)
        selected_set.add(candidate_id)
        if len(selected) == target_count:
            break

    if len(selected) != target_count:
        raise RuntimeError("Could not construct the requested bootstrap seed.")

    folds = identity_pool.set_index("candidate_id").loc[
        selected, "crossfit_fold"
    ]
    observed = set(int(value) for value in folds.tolist())
    if observed != set(range(N_FOLDS)):
        raise RuntimeError("Bootstrap did not populate all five folds.")

    return selected


class VirtualAssayOracle:
    """Truth holder. Acquisition policies never receive this object."""

    def __init__(self, truth_frame: pd.DataFrame):
        truth = truth_frame[["candidate_id", LABEL_COLUMN]].copy()
        truth["candidate_id"] = truth["candidate_id"].map(str)
        truth[LABEL_COLUMN] = truth[LABEL_COLUMN].astype(float)

        if truth["candidate_id"].duplicated().any():
            raise ValueError("Oracle candidate IDs must be unique.")
        if not np.isfinite(truth[LABEL_COLUMN].to_numpy(dtype=float)).all():
            raise ValueError("Oracle labels must contain only finite values.")

        self._labels = dict(
            zip(
                truth["candidate_id"].tolist(),
                truth[LABEL_COLUMN].tolist(),
            )
        )
        self._source_sha256 = revealed_labels_sha256(truth)

    @property
    def source_sha256(self) -> str:
        return self._source_sha256

    def contains_all(self, candidate_ids: list[str]) -> bool:
        return all(str(candidate_id) in self._labels for candidate_id in candidate_ids)

    def reveal(self, candidate_ids: list[str] | tuple[str, ...]) -> pd.DataFrame:
        requested = [str(candidate_id) for candidate_id in candidate_ids]
        if len(set(requested)) != len(requested):
            raise ValueError("Oracle reveal request contains duplicate IDs.")

        missing = [
            candidate_id
            for candidate_id in requested
            if candidate_id not in self._labels
        ]
        if missing:
            raise KeyError(f"Oracle does not contain candidate IDs: {missing}")

        return pd.DataFrame(
            {
                "candidate_id": requested,
                LABEL_COLUMN: [
                    float(self._labels[candidate_id])
                    for candidate_id in requested
                ],
            }
        )


@dataclass(frozen=True)
class PolicyView:
    round_index: int
    measurements_spent: int
    router_mode: str
    measured_ids: tuple[str, ...]
    candidates: pd.DataFrame


class AcquisitionPolicy(Protocol):
    def select(self, view: PolicyView, batch_size: int) -> list[str]:
        ...


class IdentityHashPolicy:
    """Deterministic identity-only policy used to validate the campaign loop."""

    def __init__(self, namespace: str = "NABU_PHASE2A0_POLICY"):
        self.namespace = str(namespace)

    def select(self, view: PolicyView, batch_size: int) -> list[str]:
        if batch_size < 1:
            raise ValueError("batch_size must be positive.")

        work = view.candidates[["candidate_id"]].copy()
        work["_policy_hash"] = work["candidate_id"].map(
            lambda candidate_id: fnv1a32(
                f"{self.namespace}|{candidate_id}"
            )
        )
        work = work.sort_values(
            ["_policy_hash", "candidate_id"],
            ascending=[True, True],
        )
        return work.head(batch_size)["candidate_id"].astype(str).tolist()


@dataclass(frozen=True)
class FrozenSelection:
    round_index: int
    candidate_ids: tuple[str, ...]
    selection_sha256: str
    measurements_before_reveal: int
    router_mode_before_reveal: str


def _jsonable(value):
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        return [_jsonable(item) for item in value.tolist()]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


class CampaignSimulator:
    def __init__(
        self,
        identity_pool: pd.DataFrame,
        oracle: VirtualAssayOracle,
        output_dir: str | Path | None = None,
    ):
        pool = identity_pool[
            ["candidate_id", "mutant", "mutation_set", "crossfit_fold"]
        ].copy()
        pool["candidate_id"] = pool["candidate_id"].map(str)

        if pool["candidate_id"].duplicated().any():
            raise ValueError("Campaign identity pool must have unique IDs.")
        if not oracle.contains_all(pool["candidate_id"].tolist()):
            raise ValueError(
                "Oracle does not cover every candidate in the identity pool."
            )

        self.identity_pool = pool.reset_index(drop=True)
        self.oracle = oracle
        self.output_dir = Path(output_dir) if output_dir is not None else None
        if self.output_dir is not None:
            self.output_dir.mkdir(parents=True, exist_ok=True)

        self.measured = pd.DataFrame(
            columns=[
                "candidate_id",
                "mutant",
                "mutation_set",
                "crossfit_fold",
                LABEL_COLUMN,
            ]
        )
        self.model: NabuV83Model | None = None
        self.round_index = 0
        self.bootstrap_manifest: dict | None = None
        self.round_manifests: list[dict] = []

    @property
    def measurements_spent(self) -> int:
        return int(len(self.measured))

    def _write_json(self, relative_path: str, payload: dict) -> None:
        if self.output_dir is None:
            return
        path = self.output_dir / relative_path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            json.dumps(_jsonable(payload), indent=2),
            encoding="utf-8",
        )

    def _fit(self) -> None:
        if self.measured.empty:
            raise RuntimeError("Cannot fit before any assay labels are revealed.")

        self.model = NabuV83Model().fit(
            self.measured["mutation_set"].tolist(),
            self.measured[LABEL_COLUMN].to_numpy(dtype=float),
            self.measured["candidate_id"].astype(str).tolist(),
        )

    def _diagnostics(self) -> dict:
        if self.model is None or self.model.model is None:
            raise RuntimeError("Campaign model is not fitted.")
        return {
            "router": self.model.router_decision,
            "model": model_diagnostics(self.model.model),
        }

    def bootstrap(
        self,
        seed_count: int,
        min_per_fold: int = 1,
    ) -> dict:
        if not self.measured.empty:
            raise RuntimeError("Campaign has already been bootstrapped.")

        selected_ids = bootstrap_seed_ids(
            self.identity_pool,
            target_count=seed_count,
            min_per_fold=min_per_fold,
        )
        selection_hash = selection_sha256(selected_ids)
        seed_identity = self.identity_pool[
            self.identity_pool["candidate_id"].isin(selected_ids)
        ].copy()
        seed_identity["_selection_order"] = seed_identity["candidate_id"].map(
            {candidate_id: index for index, candidate_id in enumerate(selected_ids)}
        )
        seed_identity = seed_identity.sort_values("_selection_order").drop(
            columns=["_selection_order"]
        )

        fold_counts = {
            str(fold): int(
                (seed_identity["crossfit_fold"] == fold).sum()
            )
            for fold in range(N_FOLDS)
        }

        pre_reveal = {
            "stage": "bootstrap_selection_frozen",
            "frozen_core_commit": FROZEN_CORE_COMMIT,
            "seed_count": int(len(selected_ids)),
            "candidate_ids": selected_ids,
            "selection_sha256_before_label_reveal": selection_hash,
            "crossfit_fold_counts": fold_counts,
            "truth_loaded_before_selection_freeze": False,
        }
        self._write_json("bootstrap/selection_before_reveal.json", pre_reveal)

        revealed = self.oracle.reveal(selected_ids)
        self.measured = seed_identity.merge(
            revealed,
            on="candidate_id",
            how="inner",
            validate="one_to_one",
        )
        if len(self.measured) != len(selected_ids):
            raise RuntimeError("Bootstrap reveal count mismatch.")

        self._fit()
        manifest = {
            **pre_reveal,
            "status": "BOOTSTRAP_REVEALED_AND_MODEL_FIT",
            "revealed_labels_sha256": revealed_labels_sha256(revealed),
            "measurements_after_reveal": self.measurements_spent,
            "diagnostics_after_fit": self._diagnostics(),
        }
        self.bootstrap_manifest = manifest
        self._write_json("bootstrap/manifest.json", manifest)
        self.write_campaign_manifest()
        return manifest

    def _unmeasured_identity_pool(self) -> pd.DataFrame:
        measured_ids = set(self.measured["candidate_id"].astype(str))
        return self.identity_pool[
            ~self.identity_pool["candidate_id"].isin(measured_ids)
        ].copy()

    def policy_view(self) -> PolicyView:
        if self.model is None:
            raise RuntimeError("Bootstrap the campaign before acquiring.")

        unmeasured = self._unmeasured_identity_pool()
        if unmeasured.empty:
            candidates = unmeasured.copy()
        else:
            scored = self.model.score_candidates(
                unmeasured["mutation_set"].tolist(),
                unmeasured["candidate_id"].astype(str).tolist(),
            ).drop(columns=["mutation_set"])
            candidates = unmeasured.merge(
                scored,
                on="candidate_id",
                how="left",
                validate="one_to_one",
            )

        forbidden_columns = {LABEL_COLUMN, "DMS_score", "truth", "fitness"}
        leaked = forbidden_columns.intersection(candidates.columns)
        if leaked:
            raise RuntimeError(
                f"Policy view contains forbidden truth columns: {sorted(leaked)}"
            )

        return PolicyView(
            round_index=int(self.round_index),
            measurements_spent=self.measurements_spent,
            router_mode=str(self.model.router_decision["mode"]),
            measured_ids=tuple(
                self.measured["candidate_id"].astype(str).tolist()
            ),
            candidates=candidates.copy(deep=True),
        )

    def prepare_selection(
        self,
        policy: AcquisitionPolicy,
        batch_size: int,
    ) -> FrozenSelection:
        view = self.policy_view()
        if len(view.candidates) == 0:
            raise RuntimeError("No unmeasured candidates remain.")

        actual_batch_size = min(int(batch_size), len(view.candidates))
        if actual_batch_size < 1:
            raise ValueError("batch_size must be positive.")

        selected_ids = [
            str(candidate_id)
            for candidate_id in policy.select(view, actual_batch_size)
        ]
        if len(selected_ids) != actual_batch_size:
            raise ValueError(
                "Policy must return exactly the requested batch size."
            )
        if len(set(selected_ids)) != len(selected_ids):
            raise ValueError("Policy returned duplicate candidate IDs.")

        allowed_ids = set(view.candidates["candidate_id"].astype(str))
        invalid = [
            candidate_id
            for candidate_id in selected_ids
            if candidate_id not in allowed_ids
        ]
        if invalid:
            raise ValueError(
                f"Policy selected candidates outside the unmeasured pool: {invalid}"
            )

        frozen = FrozenSelection(
            round_index=int(self.round_index),
            candidate_ids=tuple(selected_ids),
            selection_sha256=selection_sha256(selected_ids),
            measurements_before_reveal=self.measurements_spent,
            router_mode_before_reveal=str(view.router_mode),
        )

        self._write_json(
            f"rounds/round_{self.round_index:03d}_selection_before_reveal.json",
            {
                "stage": "ROUND_SELECTION_FROZEN",
                "round_index": frozen.round_index,
                "candidate_ids": list(frozen.candidate_ids),
                "selection_sha256_before_label_reveal": frozen.selection_sha256,
                "measurements_before_reveal": frozen.measurements_before_reveal,
                "router_mode_before_reveal": frozen.router_mode_before_reveal,
                "truth_loaded_before_selection_freeze": False,
            },
        )
        return frozen

    def reveal_and_update(self, frozen: FrozenSelection) -> dict:
        if frozen.round_index != self.round_index:
            raise ValueError(
                "Frozen selection round does not match current campaign round."
            )
        if selection_sha256(frozen.candidate_ids) != frozen.selection_sha256:
            raise RuntimeError("Frozen selection hash changed before reveal.")
        if self.measurements_spent != frozen.measurements_before_reveal:
            raise RuntimeError(
                "Campaign measurements changed after selection freeze."
            )

        revealed = self.oracle.reveal(frozen.candidate_ids)
        revealed_hash = revealed_labels_sha256(revealed)

        selected_identity = self.identity_pool[
            self.identity_pool["candidate_id"].isin(frozen.candidate_ids)
        ].copy()
        selected_identity["_selection_order"] = selected_identity[
            "candidate_id"
        ].map(
            {
                candidate_id: index
                for index, candidate_id in enumerate(frozen.candidate_ids)
            }
        )
        selected_identity = selected_identity.sort_values(
            "_selection_order"
        ).drop(columns=["_selection_order"])

        addition = selected_identity.merge(
            revealed,
            on="candidate_id",
            how="inner",
            validate="one_to_one",
        )
        if len(addition) != len(frozen.candidate_ids):
            raise RuntimeError("Round reveal count mismatch.")

        self.measured = pd.concat(
            [self.measured, addition],
            ignore_index=True,
        )
        if self.measured["candidate_id"].duplicated().any():
            raise RuntimeError("Campaign measured set contains duplicate IDs.")

        self._fit()
        manifest = {
            "stage": "ROUND_REVEALED_AND_REFIT",
            "round_index": int(self.round_index),
            "candidate_ids": list(frozen.candidate_ids),
            "selection_sha256_before_label_reveal": frozen.selection_sha256,
            "revealed_labels_sha256": revealed_hash,
            "measurements_before_reveal": int(
                frozen.measurements_before_reveal
            ),
            "measurements_after_reveal": self.measurements_spent,
            "router_mode_before_reveal": frozen.router_mode_before_reveal,
            "router_mode_after_refit": str(
                self.model.router_decision["mode"]
            ),
            "diagnostics_after_refit": self._diagnostics(),
        }
        self.round_manifests.append(manifest)
        self._write_json(
            f"rounds/round_{self.round_index:03d}_manifest.json",
            manifest,
        )
        self.round_index += 1
        self.write_campaign_manifest()
        return manifest

    def run_round(
        self,
        policy: AcquisitionPolicy,
        batch_size: int,
    ) -> dict:
        frozen = self.prepare_selection(policy, batch_size)
        return self.reveal_and_update(frozen)

    def write_campaign_manifest(self) -> dict:
        payload = {
            "version": "NABU_PHASE2A0_CAMPAIGN_SIMULATOR_V1",
            "frozen_core_commit": FROZEN_CORE_COMMIT,
            "oracle_source_sha256": self.oracle.source_sha256,
            "pool_size": int(len(self.identity_pool)),
            "measurements_spent": self.measurements_spent,
            "rounds_completed": int(len(self.round_manifests)),
            "bootstrap": self.bootstrap_manifest,
            "rounds": self.round_manifests,
        }
        self._write_json("CAMPAIGN_MANIFEST.json", payload)
        return payload
