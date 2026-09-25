from __future__ import annotations

import numpy as np
import pandas as pd

from baselines import _SupportPolicyBase


def percentile_rank(values: pd.Series) -> pd.Series:
    return values.rank(method="average", pct=True)


def hierarchy_rank_disagreement(frame: pd.DataFrame) -> pd.Series:
    columns = [
        "B3_RAW_PAIR",
        "B4_CROSSFIT_TRIPLET",
        "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER",
    ]
    ranks = []
    for column in columns:
        ranks.append(percentile_rank(frame[column].astype(float)))

    rank_frame = pd.concat(ranks, axis=1)
    rank_frame.columns = ["b3_rank", "b4_rank", "b5_rank"]

    disagreement = pd.concat(
        [
            (rank_frame["b3_rank"] - rank_frame["b4_rank"]).abs(),
            (rank_frame["b3_rank"] - rank_frame["b5_rank"]).abs(),
            (rank_frame["b4_rank"] - rank_frame["b5_rank"]).abs(),
        ],
        axis=1,
    ).max(axis=1)

    return disagreement.fillna(0.0).clip(lower=0.0, upper=1.0)


def state_exploration_pressure(
    support_gap: float,
    scoreability_gap: float,
    disagreement_gap: float,
) -> float:
    gaps = np.asarray(
        [support_gap, scoreability_gap, disagreement_gap],
        dtype=float,
    )
    if not np.isfinite(gaps).all():
        raise ValueError("Controller state gaps must be finite.")
    if ((gaps < 0.0) | (gaps > 1.0)).any():
        raise ValueError("Controller state gaps must be in [0, 1].")

    knowledge_completeness = float(np.prod(1.0 - gaps))
    return float(np.clip(1.0 - knowledge_completeness, 0.0, 1.0))


class EvidenceAdaptivePolicyV1(_SupportPolicyBase):
    name = "evidence_adaptive_v1"

    def __init__(self, identity_pool: pd.DataFrame):
        super().__init__(identity_pool)
        self.last_diagnostics: dict | None = None
        self.diagnostic_history: list[dict] = []

    def select(self, view, batch_size: int) -> list[str]:
        if batch_size < 1:
            raise ValueError("batch_size must be positive.")

        work = view.candidates.copy()
        if work.empty:
            raise ValueError("No candidates available for acquisition.")

        structural = self.structural_exploration(view).clip(
            lower=0.0,
            upper=1.0,
        )
        structural_rank = percentile_rank(structural)

        scoreable = work["scoreable"].astype(bool)
        scoreability_gap = float(1.0 - scoreable.mean())

        disagreement = hierarchy_rank_disagreement(work)
        disagreement_gap = float(disagreement.mean())

        support_gap = float(structural.mean())
        pressure = state_exploration_pressure(
            support_gap=support_gap,
            scoreability_gap=scoreability_gap,
            disagreement_gap=disagreement_gap,
        )

        exploitation_rank = percentile_rank(
            work["V8_3_ADAPTIVE_ROUTER"].astype(float)
        ).fillna(0.0)

        unscoreable = (~scoreable).astype(float)
        exploration_evidence = pd.concat(
            [
                structural_rank.rename("structural_rank"),
                disagreement.rename("hierarchy_disagreement"),
                unscoreable.rename("unscoreable"),
            ],
            axis=1,
        ).max(axis=1)

        acquisition_score = (
            pressure * exploration_evidence
            + (1.0 - pressure) * exploitation_rank
        )

        ranked = pd.DataFrame(
            {
                "candidate_id": work["candidate_id"].astype(str),
                "acquisition_score": acquisition_score,
                "exploration_evidence": exploration_evidence,
                "exploitation_rank": exploitation_rank,
            },
            index=work.index,
        ).sort_values(
            ["acquisition_score", "candidate_id"],
            ascending=[False, True],
        )

        selected = ranked.head(batch_size)["candidate_id"].tolist()

        self.last_diagnostics = {
            "round_index": int(view.round_index),
            "measurements_spent": int(view.measurements_spent),
            "router_mode": str(view.router_mode),
            "support_gap": support_gap,
            "scoreability_gap": scoreability_gap,
            "disagreement_gap": disagreement_gap,
            "exploration_pressure": pressure,
            "exploitation_pressure": float(1.0 - pressure),
            "candidate_count": int(len(work)),
            "scoreable_count": int(scoreable.sum()),
            "selected_ids": selected,
            "selected_mean_exploration_evidence": float(
                ranked.head(batch_size)["exploration_evidence"].mean()
            ),
            "selected_mean_exploitation_rank": float(
                ranked.head(batch_size)["exploitation_rank"].mean()
            ),
        }
        self.diagnostic_history.append(dict(self.last_diagnostics))
        return selected


def relative_exploration_pressure(
    current_gap: float,
    reference_gap: float,
) -> float:
    values = np.asarray([current_gap, reference_gap], dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("Relative controller gaps must be finite.")
    if ((values < 0.0) | (values > 1.0)).any():
        raise ValueError("Relative controller gaps must be in [0, 1].")

    denominator = float(current_gap + reference_gap)
    if denominator <= 0.0:
        return 0.0
    return float(np.clip(current_gap / denominator, 0.0, 1.0))


class RelativeEvidencePolicyV2(_SupportPolicyBase):
    name = "relative_evidence_v2"

    def __init__(self, identity_pool: pd.DataFrame):
        super().__init__(identity_pool)
        self.reference_gap: float | None = None
        self.last_diagnostics: dict | None = None
        self.diagnostic_history: list[dict] = []

    def select(self, view, batch_size: int) -> list[str]:
        if batch_size < 1:
            raise ValueError("batch_size must be positive.")

        work = view.candidates.copy()
        if work.empty:
            raise ValueError("No candidates available for acquisition.")

        structural = self.structural_exploration(view).clip(
            lower=0.0,
            upper=1.0,
        )
        structural_rank = percentile_rank(structural)

        scoreable = work["scoreable"].astype(bool)
        scoreability_gap = float(1.0 - scoreable.mean())

        disagreement = hierarchy_rank_disagreement(work)
        disagreement_gap = float(disagreement.mean())

        support_gap = float(structural.mean())
        current_gap = state_exploration_pressure(
            support_gap=support_gap,
            scoreability_gap=scoreability_gap,
            disagreement_gap=disagreement_gap,
        )

        if self.reference_gap is None:
            self.reference_gap = float(current_gap)

        pressure = relative_exploration_pressure(
            current_gap=current_gap,
            reference_gap=self.reference_gap,
        )

        exploitation_rank = percentile_rank(
            work["V8_3_ADAPTIVE_ROUTER"].astype(float)
        ).fillna(0.0)

        exploration_evidence = pd.concat(
            [
                structural_rank.rename("structural_rank"),
                disagreement.rename("hierarchy_disagreement"),
            ],
            axis=1,
        ).max(axis=1)

        acquisition_score = (
            pressure * exploration_evidence
            + (1.0 - pressure) * exploitation_rank
        )

        ranked = pd.DataFrame(
            {
                "candidate_id": work["candidate_id"].astype(str),
                "acquisition_score": acquisition_score,
                "exploration_evidence": exploration_evidence,
                "exploitation_rank": exploitation_rank,
            },
            index=work.index,
        ).sort_values(
            ["acquisition_score", "candidate_id"],
            ascending=[False, True],
        )

        selected = ranked.head(batch_size)["candidate_id"].tolist()

        self.last_diagnostics = {
            "round_index": int(view.round_index),
            "measurements_spent": int(view.measurements_spent),
            "router_mode": str(view.router_mode),
            "support_gap": support_gap,
            "scoreability_gap": scoreability_gap,
            "disagreement_gap": disagreement_gap,
            "current_gap": float(current_gap),
            "reference_gap": float(self.reference_gap),
            "exploration_pressure": pressure,
            "exploitation_pressure": float(1.0 - pressure),
            "candidate_count": int(len(work)),
            "scoreable_count": int(scoreable.sum()),
            "selected_ids": selected,
            "selected_mean_exploration_evidence": float(
                ranked.head(batch_size)["exploration_evidence"].mean()
            ),
            "selected_mean_exploitation_rank": float(
                ranked.head(batch_size)["exploitation_rank"].mean()
            ),
        }
        self.diagnostic_history.append(dict(self.last_diagnostics))
        return selected
