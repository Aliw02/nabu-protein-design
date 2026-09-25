from __future__ import annotations

from dataclasses import dataclass
import hashlib
from itertools import combinations
import json
from pathlib import Path

import numpy as np
import pandas as pd

from assembly import CandidateAssembler, MutationVocabulary
from baselines import DeterministicRandomPolicy, HistoricalFiftyFiftyPolicy
from campaign import LABEL_COLUMN, PolicyView, bootstrap_seed_ids
from nabu_protein.v83 import NabuV83Model


def selection_sha256(candidate_ids) -> str:
    payload = "\n".join(str(value) for value in candidate_ids)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


class FullLandscapeOracle:
    def __init__(self, full_frame: pd.DataFrame):
        required = {"candidate_id", "mutation_set", "DMS_score"}
        missing = required - set(full_frame.columns)
        if missing:
            raise ValueError(
                f"Oracle frame missing columns: {sorted(missing)}"
            )

        frame = full_frame[
            ["candidate_id", "mutation_set", "DMS_score"]
        ].copy()
        frame["candidate_id"] = frame["candidate_id"].astype(str)
        frame["DMS_score"] = frame["DMS_score"].astype(float)

        if frame["candidate_id"].duplicated().any():
            raise ValueError("Oracle candidate IDs must be unique.")
        if not np.isfinite(frame["DMS_score"]).all():
            raise ValueError("Oracle labels must be finite.")

        self.identity_by_id = dict(
            zip(frame["candidate_id"], frame["mutation_set"])
        )
        self.label_by_id = dict(
            zip(frame["candidate_id"], frame["DMS_score"])
        )
        self.ids = frozenset(self.label_by_id)

    def reveal(self, candidate_ids) -> pd.DataFrame:
        candidate_ids = [str(value) for value in candidate_ids]
        if len(set(candidate_ids)) != len(candidate_ids):
            raise ValueError("Oracle reveal contains duplicate candidate IDs.")
        missing = [
            candidate_id
            for candidate_id in candidate_ids
            if candidate_id not in self.label_by_id
        ]
        if missing:
            raise KeyError(f"Oracle missing candidate IDs: {missing[:5]}")

        return pd.DataFrame(
            {
                "candidate_id": candidate_ids,
                LABEL_COLUMN: [
                    float(self.label_by_id[candidate_id])
                    for candidate_id in candidate_ids
                ],
            }
        )


@dataclass(frozen=True)
class GateDiagnostics:
    ready: bool
    probe_count: int
    mature_count: int
    required_batch_size: int
    min_main_support: int
    min_pair_support: int
    mature_candidate_ids: tuple[str, ...]
    assembly_trace: tuple[dict, ...]


class EvidenceMaturityGate:
    def __init__(
        self,
        vocabulary: MutationVocabulary,
        reservoir_ids,
        assayable_ids,
        batch_size: int = 8,
        min_main_support: int = 2,
        min_pair_support: int = 2,
        beam_width: int = 256,
        probe_count: int = 64,
        target_order: int = 3,
    ):
        self.vocabulary = vocabulary
        self.reservoir_ids = frozenset(str(value) for value in reservoir_ids)
        self.assayable_ids = frozenset(str(value) for value in assayable_ids)
        self.batch_size = int(batch_size)
        self.min_main_support = int(min_main_support)
        self.min_pair_support = int(min_pair_support)
        self.beam_width = int(beam_width)
        self.probe_count = int(probe_count)
        self.target_order = int(target_order)

        if self.batch_size < 1:
            raise ValueError("batch_size must be positive.")
        if self.min_main_support < 1 or self.min_pair_support < 1:
            raise ValueError("Support thresholds must be positive.")

    def _support_mature(self, mutation_set, model) -> bool:
        main = model.model["base"]["main"]
        pair = model.model["base"]["pair"]

        for mutation in mutation_set:
            entry = main.get(mutation)
            if entry is None or int(entry["support"]) < self.min_main_support:
                return False

        for subset in combinations(mutation_set, 2):
            entry = pair.get(subset)
            if entry is None or int(entry["support"]) < self.min_pair_support:
                return False

        return True

    def probe(
        self,
        model: NabuV83Model,
        measured_mutation_sets,
        measured_ids,
    ) -> tuple[GateDiagnostics, pd.DataFrame]:
        assembler = CandidateAssembler(
            model=model,
            vocabulary=self.vocabulary,
            measured_mutation_sets=measured_mutation_sets,
        )
        proposals, trace = assembler.assemble(
            target_order=self.target_order,
            beam_width=self.beam_width,
            proposal_count=self.probe_count,
        )

        measured_ids = set(str(value) for value in measured_ids)
        mask = []
        for row in proposals.itertuples(index=False):
            candidate_id = str(row.candidate_id)
            mutation_set = tuple(row.mutation_set)
            eligible = (
                candidate_id not in measured_ids
                and candidate_id not in self.reservoir_ids
                and candidate_id in self.assayable_ids
                and bool(row.scoreable)
                and self._support_mature(mutation_set, model)
            )
            mask.append(bool(eligible))

        mature = proposals.loc[mask].copy()
        mature_ids = tuple(
            mature["candidate_id"].astype(str).tolist()
        )

        diagnostics = GateDiagnostics(
            ready=len(mature_ids) >= self.batch_size,
            probe_count=int(len(proposals)),
            mature_count=int(len(mature_ids)),
            required_batch_size=self.batch_size,
            min_main_support=self.min_main_support,
            min_pair_support=self.min_pair_support,
            mature_candidate_ids=mature_ids,
            assembly_trace=tuple(trace),
        )
        return diagnostics, mature


class DesignCampaign:
    def __init__(
        self,
        full_frame: pd.DataFrame,
        reservoir: pd.DataFrame,
        output_dir: Path | None = None,
    ):
        self.full_frame = full_frame.copy()
        self.reservoir = reservoir.copy()
        self.oracle = FullLandscapeOracle(full_frame)
        self.output_dir = None if output_dir is None else Path(output_dir)

        self.full_identity = dict(
            zip(
                self.full_frame["candidate_id"].astype(str),
                self.full_frame["mutation_set"],
            )
        )
        self.reservoir_ids = frozenset(
            self.reservoir["candidate_id"].astype(str)
        )

        self.measured = pd.DataFrame(
            columns=[
                "candidate_id",
                "mutation_set",
                LABEL_COLUMN,
                "measurement_source",
            ]
        )
        self.model: NabuV83Model | None = None
        self.round_index = 0
        self.round_manifests: list[dict] = []
        self.bootstrap_ids: tuple[str, ...] = ()

    @property
    def measurements_spent(self) -> int:
        return int(len(self.measured))

    def _write_json(self, relative_path: str, payload: dict) -> None:
        if self.output_dir is None:
            return
        path = self.output_dir / relative_path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    def _fit(self) -> None:
        self.model = NabuV83Model().fit(
            self.measured["mutation_set"].tolist(),
            self.measured[LABEL_COLUMN].to_numpy(dtype=float),
            self.measured["candidate_id"].astype(str).tolist(),
        )

    def _append_reveal(self, candidate_ids, source: str) -> None:
        candidate_ids = [str(value) for value in candidate_ids]
        already = set(self.measured["candidate_id"].astype(str))
        overlap = sorted(set(candidate_ids) & already)
        if overlap:
            raise ValueError(
                f"Attempted to measure candidates twice: {overlap[:5]}"
            )

        revealed = self.oracle.reveal(candidate_ids)
        identity = pd.DataFrame(
            {
                "candidate_id": candidate_ids,
                "mutation_set": [
                    self.full_identity[candidate_id]
                    for candidate_id in candidate_ids
                ],
            }
        )
        batch = identity.merge(
            revealed,
            on="candidate_id",
            how="inner",
            validate="one_to_one",
        )
        batch["measurement_source"] = str(source)
        self.measured = pd.concat(
            [self.measured, batch],
            ignore_index=True,
        )
        self._fit()

    def bootstrap(self, seed_count: int) -> dict:
        pool = self.reservoir[
            ["candidate_id", "crossfit_fold"]
        ].copy()
        selected = bootstrap_seed_ids(
            pool,
            target_count=int(seed_count),
            min_per_fold=1,
            namespace="NABU_PHASE2C_GB1_BOOTSTRAP",
        )
        self.bootstrap_ids = tuple(selected)

        manifest = {
            "stage": "BOOTSTRAP_SELECTION_FROZEN",
            "candidate_ids": selected,
            "selection_sha256_before_label_reveal": selection_sha256(selected),
            "truth_exposed_before_freeze": False,
        }
        self._write_json("bootstrap/selection_before_reveal.json", manifest)
        self._append_reveal(selected, source="bootstrap")
        manifest["measurements_after_reveal"] = self.measurements_spent
        manifest["router_mode_after_fit"] = str(
            self.model.router_decision["mode"]
        )
        self._write_json("bootstrap/manifest.json", manifest)
        return manifest

    def policy_view(self) -> PolicyView:
        if self.model is None:
            raise RuntimeError("Campaign is not bootstrapped.")

        measured_ids = set(self.measured["candidate_id"].astype(str))
        candidates = self.reservoir[
            ~self.reservoir["candidate_id"].astype(str).isin(measured_ids)
        ].copy()

        if len(candidates):
            scored = self.model.score_candidates(
                candidates["mutation_set"].tolist(),
                candidates["candidate_id"].astype(str).tolist(),
            ).drop(columns=["mutation_set"])
            candidates = candidates.merge(
                scored,
                on="candidate_id",
                how="left",
                validate="one_to_one",
            )

        return PolicyView(
            round_index=int(self.round_index),
            measurements_spent=self.measurements_spent,
            router_mode=str(self.model.router_decision["mode"]),
            measured_ids=tuple(
                self.measured["candidate_id"].astype(str).tolist()
            ),
            candidates=candidates,
        )

    def acquisition_selection(self, policy, batch_size: int) -> list[str]:
        view = self.policy_view()
        actual = min(int(batch_size), len(view.candidates))
        return [
            str(value)
            for value in policy.select(view, actual)
        ]

    def reveal_frozen_batch(
        self,
        candidate_ids,
        source: str,
        gate: GateDiagnostics | None = None,
    ) -> dict:
        candidate_ids = [str(value) for value in candidate_ids]
        manifest = {
            "stage": "ROUND_SELECTION_FROZEN",
            "round_index": int(self.round_index),
            "measurement_source": str(source),
            "measurements_before_reveal": self.measurements_spent,
            "candidate_ids": candidate_ids,
            "selection_sha256_before_label_reveal": selection_sha256(
                candidate_ids
            ),
            "truth_exposed_before_freeze": False,
        }
        if gate is not None:
            manifest["gate"] = {
                "ready": bool(gate.ready),
                "probe_count": int(gate.probe_count),
                "mature_count": int(gate.mature_count),
                "required_batch_size": int(gate.required_batch_size),
                "min_main_support": int(gate.min_main_support),
                "min_pair_support": int(gate.min_pair_support),
            }

        self._write_json(
            f"rounds/round_{self.round_index:03d}_selection_before_reveal.json",
            manifest,
        )
        self._append_reveal(candidate_ids, source=source)
        manifest["measurements_after_reveal"] = self.measurements_spent
        manifest["router_mode_after_fit"] = str(
            self.model.router_decision["mode"]
        )
        self.round_manifests.append(dict(manifest))
        self._write_json(
            f"rounds/round_{self.round_index:03d}_manifest.json",
            manifest,
        )
        self.round_index += 1
        return manifest


class FullLandscapeEvaluator:
    def __init__(self, full_frame: pd.DataFrame, bootstrap_count: int):
        frame = full_frame[
            ["candidate_id", "mutation_set", "DMS_score"]
        ].copy()
        self.truth_by_id = dict(
            zip(
                frame["candidate_id"].astype(str),
                frame["DMS_score"].astype(float),
            )
        )
        values = frame["DMS_score"].to_numpy(dtype=float)
        self.global_best = float(np.max(values))
        self.global_min = float(np.min(values))
        self.top1_cutoff = float(np.quantile(values, 0.99))
        self.bootstrap_count = int(bootstrap_count)

    def checkpoint(self, campaign: DesignCampaign) -> dict:
        measured = campaign.measured.copy()
        values = measured[LABEL_COLUMN].to_numpy(dtype=float)
        best = float(np.max(values))
        denominator = self.global_best - self.global_min
        normalized_best = (
            1.0
            if denominator <= 0
            else float((best - self.global_min) / denominator)
        )
        top_mask = values >= self.top1_cutoff

        top_sources = (
            measured.loc[top_mask, "measurement_source"]
            .value_counts()
            .to_dict()
        )

        return {
            "measurements_spent": campaign.measurements_spent,
            "best_true_fitness": best,
            "normalized_best_fitness": normalized_best,
            "best_so_far_regret": float(self.global_best - best),
            "cumulative_top1_hits": int(np.sum(top_mask)),
            "top1_hits_by_source": {
                str(key): int(value)
                for key, value in top_sources.items()
            },
            "assembly_measurements": int(
                (measured["measurement_source"] == "assembly").sum()
            ),
            "assembly_top1_hits": int(
                (
                    top_mask
                    & (
                        measured["measurement_source"].to_numpy()
                        == "assembly"
                    )
                ).sum()
            ),
            "router_mode": str(campaign.model.router_decision["mode"]),
            "visible_b3_oof_spearman": float(
                campaign.model.router_decision["B3_oof_spearman"]
            ),
        }

    def summarize(self, curve: list[dict]) -> dict:
        frame = pd.DataFrame(curve).sort_values("measurements_spent")
        post = frame[
            frame["measurements_spent"] >= self.bootstrap_count
        ].copy()

        if len(post) < 2:
            auc = float(post["normalized_best_fitness"].iloc[0])
        else:
            x = post["measurements_spent"].to_numpy(dtype=float)
            y = post["normalized_best_fitness"].to_numpy(dtype=float)
            span = float(x[-1] - x[0])
            auc = (
                float(y[-1])
                if span <= 0
                else float(np.trapezoid(y, x) / span)
            )

        top_rows = frame[frame["cumulative_top1_hits"] > 0]
        first_top1 = (
            None
            if top_rows.empty
            else int(top_rows.iloc[0]["measurements_spent"])
        )

        return {
            "global_best_true_fitness": self.global_best,
            "top1_cutoff": self.top1_cutoff,
            "first_top1_measurement": first_top1,
            "post_bootstrap_normalized_discovery_auc": auc,
            "final": frame.iloc[-1].to_dict(),
        }


def make_policy(name: str, reservoir: pd.DataFrame):
    if name == "random":
        return DeterministicRandomPolicy(seed=161)
    if name == "historical_50_50":
        return HistoricalFiftyFiftyPolicy(reservoir)
    raise ValueError(f"Unknown acquisition policy: {name}")
