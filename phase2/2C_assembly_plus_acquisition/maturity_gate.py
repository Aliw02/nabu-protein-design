from __future__ import annotations

from itertools import combinations

import pandas as pd

from assembly import CandidateAssembler


def proposal_is_mature(
    mutation_set: tuple[str, ...],
    model,
    min_support: int = 2,
) -> bool:
    main_memory = model.model["base"]["main"]
    pair_memory = model.model["base"]["pair"]

    for mutation in mutation_set:
        entry = main_memory.get(mutation)
        if entry is None:
            return False
        if int(entry.get("support", 0)) < min_support:
            return False

    for pair in combinations(mutation_set, 2):
        entry = pair_memory.get(pair)
        if entry is None:
            return False
        if int(entry.get("support", 0)) < min_support:
            return False

    return True


def mature_assembly_frontier(
    model,
    vocabulary,
    codec,
    measured_frame: pd.DataFrame,
    reservoir_ids: set[str],
    assay_universe_ids: set[str],
    target_order: int = 3,
    beam_width: int = 256,
    proposal_frontier: int = 1024,
    min_support: int = 2,
) -> tuple[pd.DataFrame, dict]:
    assembler = CandidateAssembler(
        model=model,
        vocabulary=vocabulary,
        measured_mutation_sets=measured_frame[
            "mutation_set"
        ].tolist(),
    )
    proposals, trace = assembler.assemble(
        target_order=target_order,
        beam_width=beam_width,
        proposal_count=proposal_frontier,
    )

    measured_ids = set(
        measured_frame["candidate_id"].astype(str)
    )

    work = proposals.copy()
    work["genotype_id"] = work["mutation_set"].map(
        codec.decode
    )
    work["in_assay_universe"] = work[
        "genotype_id"
    ].isin(assay_universe_ids)
    work["outside_reservoir"] = ~work[
        "genotype_id"
    ].isin(reservoir_ids)
    work["unmeasured"] = ~work[
        "genotype_id"
    ].isin(measured_ids)
    work["support_mature"] = work[
        "mutation_set"
    ].map(
        lambda mutation_set: proposal_is_mature(
            mutation_set,
            model,
            min_support=min_support,
        )
    )

    mature_mask = (
        work["in_assay_universe"]
        & work["outside_reservoir"]
        & work["unmeasured"]
        & work["support_mature"]
        & work["scoreable"].astype(bool)
    )
    mature = work[mature_mask].copy()

    diagnostics = {
        "target_order": int(target_order),
        "beam_width": int(beam_width),
        "proposal_frontier": int(proposal_frontier),
        "min_support": int(min_support),
        "raw_proposal_count": int(len(work)),
        "assayable_count": int(work["in_assay_universe"].sum()),
        "outside_reservoir_count": int(
            work["outside_reservoir"].sum()
        ),
        "unmeasured_count": int(work["unmeasured"].sum()),
        "support_mature_count": int(
            work["support_mature"].sum()
        ),
        "mature_count": int(len(mature)),
        "assembly_trace": trace,
    }
    return mature.reset_index(drop=True), diagnostics


class EvidenceMaturityGateV1:
    name = "evidence_maturity_gate_v1"

    def __init__(self, batch_size: int = 8):
        if batch_size < 1:
            raise ValueError(
                "batch_size must be positive."
            )
        self.batch_size = int(batch_size)
        self.latched_open = False
        self.opened_at_measurement: int | None = None
        self.history: list[dict] = []

    def update(
        self,
        measurements_spent: int,
        mature_count: int,
        diagnostics: dict,
    ) -> bool:
        ready_now = int(mature_count) >= self.batch_size

        if ready_now and not self.latched_open:
            self.latched_open = True
            self.opened_at_measurement = int(
                measurements_spent
            )

        row = {
            "measurements_spent": int(
                measurements_spent
            ),
            "ready_now": bool(ready_now),
            "latched_open": bool(
                self.latched_open
            ),
            "opened_at_measurement": (
                None
                if self.opened_at_measurement is None
                else int(
                    self.opened_at_measurement
                )
            ),
            **diagnostics,
        }
        self.history.append(row)
        return bool(self.latched_open)


class EvidenceMaturityGateV3:
    """Per-round dynamic maturity gate for Phase 2C V3."""

    name = "evidence_maturity_gate_v3"

    def __init__(self, batch_size: int = 8):
        if batch_size < 1:
            raise ValueError("batch_size must be positive.")
        self.batch_size = int(batch_size)
        self.first_ready_at_measurement: int | None = None
        self.ready_rounds = 0
        self.not_ready_rounds = 0
        self.fallback_rounds_after_first_ready = 0
        self.transition_count = 0
        self._previous_ready: bool | None = None
        self.history: list[dict] = []

    def update(
        self,
        measurements_spent: int,
        mature_count: int,
        diagnostics: dict,
    ) -> bool:
        ready_now = int(mature_count) >= self.batch_size

        if ready_now:
            self.ready_rounds += 1
            if self.first_ready_at_measurement is None:
                self.first_ready_at_measurement = int(
                    measurements_spent
                )
        else:
            self.not_ready_rounds += 1
            if self.first_ready_at_measurement is not None:
                self.fallback_rounds_after_first_ready += 1

        if (
            self._previous_ready is not None
            and ready_now != self._previous_ready
        ):
            self.transition_count += 1
        self._previous_ready = bool(ready_now)

        row = {
            "measurements_spent": int(measurements_spent),
            "ready_now": bool(ready_now),
            "first_ready_at_measurement": (
                None
                if self.first_ready_at_measurement is None
                else int(self.first_ready_at_measurement)
            ),
            "ready_rounds": int(self.ready_rounds),
            "not_ready_rounds": int(self.not_ready_rounds),
            "fallback_rounds_after_first_ready": int(
                self.fallback_rounds_after_first_ready
            ),
            "transition_count": int(self.transition_count),
            **diagnostics,
        }
        self.history.append(row)
        return bool(ready_now)
