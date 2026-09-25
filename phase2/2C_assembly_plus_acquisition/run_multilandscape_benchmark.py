from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
PHASE2 = HERE.parent
PHASE2A0 = (
    PHASE2
    / "2A_closed_loop_acquisition"
    / "2A0_campaign_simulator"
)
PHASE2A1 = (
    PHASE2
    / "2A_closed_loop_acquisition"
    / "2A1_baselines"
)
PHASE2B = PHASE2 / "2B_candidate_assembly"

for path in (PHASE2A0, PHASE2A1, PHASE2B):
    sys.path.insert(0, str(path))

from campaign import (  # noqa: E402
    LABEL_COLUMN,
    PolicyView,
    VirtualAssayOracle,
)
from baselines import (  # noqa: E402
    DeterministicRandomPolicy,
    HistoricalFiftyFiftyPolicy,
)
from assembly import (  # noqa: E402
    MutationVocabulary,
    ReferenceProtein,
)
from maturity_gate import (  # noqa: E402
    EvidenceMaturityGateV1,
    mature_assembly_frontier,
)
from four_site_codec import (  # noqa: E402
    bootstrap_ids,
    deterministic_identity_subpool,
    load_four_site_landscape,
    mutation_vocabulary_from_identity,
    sha256_ids,
)
from nabu_protein.higher_order import model_diagnostics  # noqa: E402
from nabu_protein.v83 import NabuV83Model  # noqa: E402


LANDSCAPE_SPECS = {
    "GB1": {
        "variant_column": "variant",
        "fitness_column": "fitness",
        "reference_genotype": "VDGV",
    },
    "TRPB": {
        "variant_column": "variant",
        "fitness_column": "fitness",
        "reference_genotype": "VFVS",
    },
    "PHOQ": {
        "variant_column": "Variants",
        "fitness_column": "Fitness",
        "reference_genotype": "AVST",
    },
}

CONDITIONS = (
    "random_only",
    "acquisition_only",
    "random_then_gated_assembly",
    "acquisition_then_gated_assembly",
)


def _json_default(value):
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.bool_):
        return bool(value)
    raise TypeError(
        f"Object of type {type(value).__name__} is not JSON serializable."
    )


def _file_sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _truth_frame(universe: pd.DataFrame) -> pd.DataFrame:
    return universe[
        ["candidate_id", "DMS_score"]
    ].rename(
        columns={"DMS_score": LABEL_COLUMN}
    )


def _fit_model(measured: pd.DataFrame) -> NabuV83Model:
    return NabuV83Model().fit(
        measured["mutation_set"].tolist(),
        measured[LABEL_COLUMN].to_numpy(dtype=float),
        measured["candidate_id"].astype(str).tolist(),
    )


def _score_identity_frame(
    model: NabuV83Model,
    frame: pd.DataFrame,
) -> pd.DataFrame:
    if frame.empty:
        return frame.copy()

    return model.score_candidates(
        frame["mutation_set"].tolist(),
        frame["candidate_id"].astype(str).tolist(),
    )


def _policy_batch(
    policy,
    model: NabuV83Model,
    measured: pd.DataFrame,
    candidates: pd.DataFrame,
    round_index: int,
    batch_size: int,
) -> list[str]:
    if candidates.empty:
        raise RuntimeError("No acquisition candidates remain.")

    scored = _score_identity_frame(
        model,
        candidates[
            ["candidate_id", "mutation_set"]
        ].copy(),
    )
    view = PolicyView(
        round_index=int(round_index),
        measurements_spent=int(len(measured)),
        router_mode=str(
            model.router_decision["mode"]
        ),
        measured_ids=tuple(
            measured["candidate_id"].astype(str)
        ),
        candidates=scored,
    )
    return [
        str(candidate_id)
        for candidate_id in policy.select(
            view,
            min(batch_size, len(scored)),
        )
    ]


def _reveal_into_measured(
    oracle: VirtualAssayOracle,
    identity_lookup: pd.DataFrame,
    measured: pd.DataFrame,
    selected_ids: list[str],
) -> pd.DataFrame:
    revealed = oracle.reveal(selected_ids)
    identity = identity_lookup[
        identity_lookup["candidate_id"].isin(selected_ids)
    ][
        ["candidate_id", "mutation_set", "crossfit_fold"]
    ].copy()

    order = {
        candidate_id: index
        for index, candidate_id in enumerate(selected_ids)
    }
    identity["_order"] = identity["candidate_id"].map(order)
    identity = identity.sort_values("_order").drop(
        columns=["_order"]
    )

    new_rows = identity.merge(
        revealed,
        on="candidate_id",
        how="inner",
        validate="one_to_one",
    )
    if len(new_rows) != len(selected_ids):
        raise RuntimeError(
            "Reveal did not map one-to-one onto selected identities."
        )

    updated = pd.concat(
        [measured, new_rows],
        ignore_index=True,
    )
    if updated["candidate_id"].duplicated().any():
        raise RuntimeError(
            "Campaign attempted to measure a candidate twice."
        )
    return updated


class FullUniverseEvaluator:
    def __init__(
        self,
        universe: pd.DataFrame,
        bootstrap_count: int,
    ):
        self.truth_by_id = dict(
            zip(
                universe["candidate_id"].astype(str),
                universe["DMS_score"].astype(float),
            )
        )
        self.mutation_by_id = dict(
            zip(
                universe["candidate_id"].astype(str),
                universe["mutation_set"],
            )
        )
        values = universe[
            "DMS_score"
        ].to_numpy(dtype=float)
        self.global_best = float(np.max(values))
        self.global_min = float(np.min(values))
        self.top1_cutoff = float(
            np.quantile(values, 0.99)
        )
        self.bootstrap_count = int(
            bootstrap_count
        )

    def checkpoint(
        self,
        measured: pd.DataFrame,
        model: NabuV83Model,
        reservoir_unmeasured: pd.DataFrame,
        assembly_selected_total: int,
        assembly_outside_reservoir_total: int,
    ) -> dict:
        measured_ids = measured[
            "candidate_id"
        ].astype(str).tolist()
        values = np.asarray(
            [
                self.truth_by_id[candidate_id]
                for candidate_id in measured_ids
            ],
            dtype=float,
        )
        best = float(np.max(values))
        top1_hits = int(
            np.sum(values >= self.top1_cutoff)
        )

        denominator = (
            self.global_best - self.global_min
        )
        normalized_best = (
            1.0
            if denominator <= 0.0
            else float(
                (best - self.global_min)
                / denominator
            )
        )

        if reservoir_unmeasured.empty:
            scoreable_coverage = 1.0
        else:
            scored = _score_identity_frame(
                model,
                reservoir_unmeasured[
                    ["candidate_id", "mutation_set"]
                ],
            )
            scoreable_coverage = float(
                scored["scoreable"].astype(bool).mean()
            )

        return {
            "measurements_spent": int(len(measured)),
            "best_true_fitness": best,
            "normalized_best_fitness": normalized_best,
            "best_so_far_regret": float(
                self.global_best - best
            ),
            "cumulative_top1_hits": top1_hits,
            "scoreable_coverage_reservoir": scoreable_coverage,
            "router_mode": str(
                model.router_decision["mode"]
            ),
            "visible_b3_oof_spearman": float(
                model.router_decision[
                    "B3_oof_spearman"
                ]
            ),
            "assembly_selected_total": int(
                assembly_selected_total
            ),
            "assembly_outside_reservoir_total": int(
                assembly_outside_reservoir_total
            ),
        }

    def summarize(
        self,
        curve: list[dict],
        bootstrap_ids_set: set[str],
        final_measured_ids: list[str],
        assembly_selected_ids: list[str],
    ) -> dict:
        frame = pd.DataFrame(curve).sort_values(
            "measurements_spent"
        )
        post = frame[
            frame["measurements_spent"]
            >= self.bootstrap_count
        ].copy()

        if len(post) < 2:
            auc = float(
                post[
                    "normalized_best_fitness"
                ].iloc[0]
            )
        else:
            x = post[
                "measurements_spent"
            ].to_numpy(dtype=float)
            y = post[
                "normalized_best_fitness"
            ].to_numpy(dtype=float)
            span = float(x[-1] - x[0])
            auc = (
                float(y[-1])
                if span <= 0.0
                else float(
                    np.trapezoid(y, x) / span
                )
            )

        top1_rows = frame[
            frame["cumulative_top1_hits"] > 0
        ]
        first_top1 = (
            None
            if top1_rows.empty
            else int(
                top1_rows.iloc[0][
                    "measurements_spent"
                ]
            )
        )

        acquired_ids = [
            candidate_id
            for candidate_id in final_measured_ids
            if candidate_id not in bootstrap_ids_set
        ]
        acquired_top1 = int(
            sum(
                self.truth_by_id[candidate_id]
                >= self.top1_cutoff
                for candidate_id in acquired_ids
            )
        )
        assembled_top1 = int(
            sum(
                self.truth_by_id[candidate_id]
                >= self.top1_cutoff
                for candidate_id in assembly_selected_ids
            )
        )

        final = frame.iloc[-1].to_dict()
        return {
            "global_best_true_fitness": self.global_best,
            "top1_cutoff": self.top1_cutoff,
            "post_bootstrap_normalized_discovery_auc": auc,
            "first_top1_measurement": first_top1,
            "acquired_top1_hits_after_bootstrap": acquired_top1,
            "assembled_top1_hits": assembled_top1,
            "final": final,
        }


def _batch_truth_diagnostic(
    oracle: VirtualAssayOracle,
    candidate_ids: list[str],
    top1_cutoff: float,
) -> dict:
    revealed = oracle.reveal(candidate_ids)
    values = revealed[
        LABEL_COLUMN
    ].to_numpy(dtype=float)
    return {
        "count": int(len(values)),
        "mean_true_fitness": float(
            np.mean(values)
        ),
        "best_true_fitness": float(
            np.max(values)
        ),
        "top1_hits": int(
            np.sum(values >= top1_cutoff)
        ),
    }


def run_condition(
    landscape_name: str,
    universe: pd.DataFrame,
    codec,
    reservoir: pd.DataFrame,
    seed_ids: list[str],
    condition: str,
    output_dir: Path | None,
    target_count: int,
    batch_size: int,
    beam_width: int,
    proposal_frontier: int,
    min_support: int,
) -> dict:
    if condition not in CONDITIONS:
        raise ValueError(
            f"Unknown condition: {condition}"
        )

    universe_identity = universe[
        ["candidate_id", "mutation_set", "crossfit_fold"]
    ].copy()
    truth = _truth_frame(universe)
    oracle = VirtualAssayOracle(truth)

    seed_identity = universe_identity[
        universe_identity["candidate_id"].isin(seed_ids)
    ].copy()
    order = {
        candidate_id: index
        for index, candidate_id in enumerate(seed_ids)
    }
    seed_identity["_order"] = seed_identity[
        "candidate_id"
    ].map(order)
    seed_identity = (
        seed_identity.sort_values("_order")
        .drop(columns=["_order"])
    )
    measured = _reveal_into_measured(
        oracle=oracle,
        identity_lookup=universe_identity,
        measured=pd.DataFrame(
            columns=[
                "candidate_id",
                "mutation_set",
                "crossfit_fold",
                LABEL_COLUMN,
            ]
        ),
        selected_ids=seed_ids,
    )
    model = _fit_model(measured)

    reference = ReferenceProtein(
        name=f"{landscape_name}_four_site_reference",
        sequence=codec.reference_genotype,
    )
    vocabulary = MutationVocabulary.from_tokens(
        reference,
        mutation_vocabulary_from_identity(
            universe,
            codec,
        ),
    )

    reservoir_ids = set(
        reservoir["candidate_id"].astype(str)
    )
    universe_ids = set(
        universe["candidate_id"].astype(str)
    )

    random_policy = DeterministicRandomPolicy(
        seed=161
    )
    acquisition_policy = HistoricalFiftyFiftyPolicy(
        universe[
            ["candidate_id", "mutation_set"]
        ]
    )
    gate = EvidenceMaturityGateV1(
        batch_size=batch_size
    )
    evaluator = FullUniverseEvaluator(
        universe=universe,
        bootstrap_count=len(seed_ids),
    )

    curve = []
    round_rows = []
    selection_hashes = []
    assembly_selected_ids = []
    assembly_active_rounds = 0
    state_matched_ablation = None

    def reservoir_unmeasured():
        measured_ids = set(
            measured["candidate_id"].astype(str)
        )
        return reservoir[
            ~reservoir["candidate_id"].isin(
                measured_ids
            )
        ].copy()

    curve.append(
        evaluator.checkpoint(
            measured=measured,
            model=model,
            reservoir_unmeasured=reservoir_unmeasured(),
            assembly_selected_total=0,
            assembly_outside_reservoir_total=0,
        )
    )

    round_index = 0
    while len(measured) < target_count:
        remaining = target_count - len(measured)
        actual_batch = min(
            batch_size,
            remaining,
        )

        current_reservoir = reservoir_unmeasured()
        mature, maturity_diagnostics = (
            mature_assembly_frontier(
                model=model,
                vocabulary=vocabulary,
                codec=codec,
                measured_frame=measured,
                reservoir_ids=reservoir_ids,
                assay_universe_ids=universe_ids,
                target_order=3,
                beam_width=beam_width,
                proposal_frontier=proposal_frontier,
                min_support=min_support,
            )
        )
        gate_open = gate.update(
            measurements_spent=len(measured),
            mature_count=len(mature),
            diagnostics=maturity_diagnostics,
        )

        use_assembly = (
            condition
            in {
                "random_then_gated_assembly",
                "acquisition_then_gated_assembly",
            }
            and gate_open
        )

        selection_source = None
        if use_assembly:
            if len(mature) < actual_batch:
                raise RuntimeError(
                    f"{landscape_name}/{condition}: "
                    "latched assembly gate cannot supply "
                    f"batch of {actual_batch}; mature={len(mature)}."
                )

            selected_ids = (
                mature.head(actual_batch)[
                    "genotype_id"
                ]
                .astype(str)
                .tolist()
            )
            selection_source = "mature_assembly"
            assembly_active_rounds += 1

            if (
                condition
                == "acquisition_then_gated_assembly"
                and state_matched_ablation is None
            ):
                alternate_ids = _policy_batch(
                    policy=HistoricalFiftyFiftyPolicy(
                        universe[
                            [
                                "candidate_id",
                                "mutation_set",
                            ]
                        ]
                    ),
                    model=model,
                    measured=measured,
                    candidates=current_reservoir,
                    round_index=round_index,
                    batch_size=actual_batch,
                )
                state_matched_ablation = {
                    "measurements_spent": int(
                        len(measured)
                    ),
                    "assembly_ids": selected_ids,
                    "assembly_sha256": sha256_ids(
                        selected_ids
                    ),
                    "acquisition_ids": alternate_ids,
                    "acquisition_sha256": sha256_ids(
                        alternate_ids
                    ),
                    "truth_exposed_before_both_freezes": False,
                }
        else:
            if condition in {
                "random_only",
                "random_then_gated_assembly",
            }:
                policy = random_policy
                selection_source = "reservoir_random"
            else:
                policy = acquisition_policy
                selection_source = "reservoir_historical_50_50"

            selected_ids = _policy_batch(
                policy=policy,
                model=model,
                measured=measured,
                candidates=current_reservoir,
                round_index=round_index,
                batch_size=actual_batch,
            )

        if len(selected_ids) != actual_batch:
            raise RuntimeError(
                "Selection did not return the requested batch."
            )
        if len(set(selected_ids)) != len(selected_ids):
            raise RuntimeError(
                "Selection contains duplicate candidate IDs."
            )
        measured_id_set = set(
            measured["candidate_id"].astype(str)
        )
        if measured_id_set.intersection(
            selected_ids
        ):
            raise RuntimeError(
                "Selection repeats a measured candidate."
            )

        selection_hash = sha256_ids(
            selected_ids
        )
        selection_hashes.append(
            selection_hash
        )

        pre_reveal = {
            "round_index": int(round_index),
            "measurements_before_reveal": int(
                len(measured)
            ),
            "selection_source": selection_source,
            "candidate_ids": selected_ids,
            "selection_sha256_before_reveal": selection_hash,
            "gate_open": bool(gate_open),
            "gate_opened_at_measurement": (
                gate.opened_at_measurement
            ),
            "maturity": maturity_diagnostics,
            "truth_exposed_to_selection_before_freeze": False,
        }

        if output_dir is not None:
            round_dir = (
                output_dir
                / "rounds"
            )
            round_dir.mkdir(
                parents=True,
                exist_ok=True,
            )
            (
                round_dir
                / f"round_{round_index:03d}_selection_before_reveal.json"
            ).write_text(
                json.dumps(
                    pre_reveal,
                    indent=2,
                    default=_json_default,
                ),
                encoding="utf-8",
            )

        if state_matched_ablation is not None:
            if (
                "assembly_truth" not in state_matched_ablation
                and state_matched_ablation[
                    "measurements_spent"
                ]
                == len(measured)
            ):
                state_matched_ablation[
                    "assembly_truth"
                ] = _batch_truth_diagnostic(
                    oracle,
                    state_matched_ablation[
                        "assembly_ids"
                    ],
                    evaluator.top1_cutoff,
                )
                state_matched_ablation[
                    "acquisition_truth"
                ] = _batch_truth_diagnostic(
                    oracle,
                    state_matched_ablation[
                        "acquisition_ids"
                    ],
                    evaluator.top1_cutoff,
                )

        measured = _reveal_into_measured(
            oracle=oracle,
            identity_lookup=universe_identity,
            measured=measured,
            selected_ids=selected_ids,
        )
        model = _fit_model(measured)

        if selection_source == "mature_assembly":
            assembly_selected_ids.extend(
                selected_ids
            )

        diagnostics = model_diagnostics(
            model.model
        )
        round_rows.append(
            {
                "round_index": int(round_index),
                "measurements_after_reveal": int(
                    len(measured)
                ),
                "selection_source": selection_source,
                "selection_sha256": selection_hash,
                "gate_open": bool(
                    gate.latched_open
                ),
                "gate_opened_at_measurement": (
                    gate.opened_at_measurement
                ),
                "mature_count_before_reveal": int(
                    maturity_diagnostics[
                        "mature_count"
                    ]
                ),
                "pair_entries_after_fit": int(
                    diagnostics["pair"]["entries"]
                ),
                "triplet_entries_after_fit": int(
                    diagnostics[
                        "triplet_crossfit"
                    ]["entries"]
                ),
                "quartet_entries_after_fit": int(
                    diagnostics[
                        "quartet_crossfit"
                    ]["entries"]
                ),
                "router_mode_after_fit": str(
                    model.router_decision["mode"]
                ),
            }
        )

        curve.append(
            evaluator.checkpoint(
                measured=measured,
                model=model,
                reservoir_unmeasured=reservoir_unmeasured(),
                assembly_selected_total=len(
                    assembly_selected_ids
                ),
                assembly_outside_reservoir_total=sum(
                    candidate_id not in reservoir_ids
                    for candidate_id in assembly_selected_ids
                ),
            )
        )
        round_index += 1

    transcript_hash = hashlib.sha256(
        "\n".join(
            selection_hashes
        ).encode("utf-8")
    ).hexdigest()

    final_ids = measured[
        "candidate_id"
    ].astype(str).tolist()
    summary = evaluator.summarize(
        curve=curve,
        bootstrap_ids_set=set(seed_ids),
        final_measured_ids=final_ids,
        assembly_selected_ids=assembly_selected_ids,
    )
    summary.update(
        {
            "version": "NABU_PHASE2C_CONDITION_RESULT_V1",
            "landscape": landscape_name,
            "condition": condition,
            "bootstrap_count": int(
                len(seed_ids)
            ),
            "target_count": int(
                target_count
            ),
            "batch_size": int(
                batch_size
            ),
            "rounds_completed": int(
                round_index
            ),
            "gate_opened_at_measurement": (
                gate.opened_at_measurement
            ),
            "assembly_active_rounds": int(
                assembly_active_rounds
            ),
            "assembly_selected_count": int(
                len(assembly_selected_ids)
            ),
            "assembly_outside_reservoir_count": int(
                sum(
                    candidate_id not in reservoir_ids
                    for candidate_id in assembly_selected_ids
                )
            ),
            "selection_transcript_sha256": transcript_hash,
            "state_matched_first_gate_ablation": state_matched_ablation,
        }
    )

    if output_dir is not None:
        output_dir.mkdir(
            parents=True,
            exist_ok=True,
        )
        pd.DataFrame(curve).to_csv(
            output_dir / "DISCOVERY_CURVE.csv",
            index=False,
        )
        pd.DataFrame(round_rows).to_csv(
            output_dir / "ROUND_SUMMARY.csv",
            index=False,
        )
        pd.DataFrame(
            gate.history
        ).to_csv(
            output_dir / "MATURITY_TRAJECTORY.csv",
            index=False,
        )
        (
            output_dir
            / "CONDITION_SUMMARY.json"
        ).write_text(
            json.dumps(
                summary,
                indent=2,
                default=_json_default,
            ),
            encoding="utf-8",
        )

    return summary


def _landscape_win(
    acquisition: dict,
    combined: dict,
) -> bool:
    return bool(
        combined[
            "post_bootstrap_normalized_discovery_auc"
        ]
        > acquisition[
            "post_bootstrap_normalized_discovery_auc"
        ]
        and combined["final"][
            "best_so_far_regret"
        ]
        <= acquisition["final"][
            "best_so_far_regret"
        ]
        and combined["final"][
            "cumulative_top1_hits"
        ]
        >= acquisition["final"][
            "cumulative_top1_hits"
        ]
        and combined[
            "assembly_active_rounds"
        ]
        > 0
    )


def _joint_degradation(
    acquisition: dict,
    combined: dict,
) -> bool:
    return bool(
        combined[
            "post_bootstrap_normalized_discovery_auc"
        ]
        < acquisition[
            "post_bootstrap_normalized_discovery_auc"
        ]
        and combined["final"][
            "best_so_far_regret"
        ]
        > acquisition["final"][
            "best_so_far_regret"
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run the preregistered NABU Phase-2C "
            "multilandscape development benchmark."
        )
    )
    parser.add_argument("--gb1", required=True)
    parser.add_argument("--trpb", required=True)
    parser.add_argument("--phoq", required=True)
    parser.add_argument(
        "--out",
        default=(
            "phase2/2C_assembly_plus_acquisition/"
            "results/multilandscape_v1"
        ),
    )
    parser.add_argument(
        "--reservoir-size",
        type=int,
        default=2048,
    )
    parser.add_argument(
        "--bootstrap-count",
        type=int,
        default=64,
    )
    parser.add_argument(
        "--target-count",
        type=int,
        default=256,
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=8,
    )
    parser.add_argument(
        "--beam-width",
        type=int,
        default=256,
    )
    parser.add_argument(
        "--proposal-frontier",
        type=int,
        default=1024,
    )
    parser.add_argument(
        "--min-support",
        type=int,
        default=2,
    )
    args = parser.parse_args()

    if not (
        0
        < args.bootstrap_count
        < args.target_count
        <= args.reservoir_size
    ):
        raise SystemExit(
            "Require 0 < bootstrap-count < target-count "
            "<= reservoir-size."
        )
    if args.batch_size < 1:
        raise SystemExit(
            "--batch-size must be positive."
        )

    input_paths = {
        "GB1": args.gb1,
        "TRPB": args.trpb,
        "PHOQ": args.phoq,
    }

    out = Path(args.out)
    out.mkdir(
        parents=True,
        exist_ok=True,
    )

    all_summaries = []
    landscape_records = []
    replay_checks = {}

    for landscape_name in (
        "GB1",
        "TRPB",
        "PHOQ",
    ):
        spec = LANDSCAPE_SPECS[
            landscape_name
        ]
        universe, codec = (
            load_four_site_landscape(
                path=input_paths[
                    landscape_name
                ],
                variant_column=spec[
                    "variant_column"
                ],
                fitness_column=spec[
                    "fitness_column"
                ],
                reference_genotype=spec[
                    "reference_genotype"
                ],
            )
        )
        reservoir = (
            deterministic_identity_subpool(
                universe,
                pool_size=args.reservoir_size,
            )
        )
        seeds = bootstrap_ids(
            reservoir,
            seed_count=args.bootstrap_count,
        )

        source_record = {
            "landscape": landscape_name,
            "input": str(
                input_paths[landscape_name]
            ),
            "source_sha256": _file_sha256(
                input_paths[landscape_name]
            ),
            "reference_genotype": spec[
                "reference_genotype"
            ],
            "assay_universe_size": int(
                len(universe)
            ),
            "reservoir_size": int(
                len(reservoir)
            ),
            "reservoir_ids_sha256": sha256_ids(
                reservoir[
                    "candidate_id"
                ].astype(str).tolist()
            ),
            "bootstrap_ids_sha256": sha256_ids(
                seeds
            ),
        }
        landscape_records.append(
            source_record
        )

        summaries = {}
        for condition in CONDITIONS:
            print(
                f"[2C] START {landscape_name} "
                f"{condition}"
            )
            condition_dir = (
                out
                / landscape_name
                / condition
            )
            summary = run_condition(
                landscape_name=landscape_name,
                universe=universe,
                codec=codec,
                reservoir=reservoir,
                seed_ids=seeds,
                condition=condition,
                output_dir=condition_dir,
                target_count=args.target_count,
                batch_size=args.batch_size,
                beam_width=args.beam_width,
                proposal_frontier=args.proposal_frontier,
                min_support=args.min_support,
            )
            summaries[condition] = summary
            all_summaries.append(
                summary
            )
            print(
                f"[2C] DONE {landscape_name} "
                f"{condition} "
                f"best={summary['final']['best_true_fitness']:.6f} "
                f"top1={int(summary['final']['cumulative_top1_hits'])} "
                f"auc={summary['post_bootstrap_normalized_discovery_auc']:.6f} "
                f"gate={summary['gate_opened_at_measurement']} "
                f"assembly={summary['assembly_selected_count']}"
            )

        replay = run_condition(
            landscape_name=landscape_name,
            universe=universe,
            codec=codec,
            reservoir=reservoir,
            seed_ids=seeds,
            condition="acquisition_then_gated_assembly",
            output_dir=None,
            target_count=args.target_count,
            batch_size=args.batch_size,
            beam_width=args.beam_width,
            proposal_frontier=args.proposal_frontier,
            min_support=args.min_support,
        )
        replay_ok = (
            replay[
                "selection_transcript_sha256"
            ]
            == summaries[
                "acquisition_then_gated_assembly"
            ][
                "selection_transcript_sha256"
            ]
        )
        if not replay_ok:
            raise RuntimeError(
                f"Deterministic replay failed for {landscape_name}."
            )
        replay_checks[
            landscape_name
        ] = True

    rows = []
    by_landscape = {}
    for summary in all_summaries:
        row = {
            "landscape": summary[
                "landscape"
            ],
            "condition": summary[
                "condition"
            ],
            "discovery_auc": summary[
                "post_bootstrap_normalized_discovery_auc"
            ],
            "final_best_true_fitness": summary[
                "final"
            ]["best_true_fitness"],
            "final_regret": summary[
                "final"
            ]["best_so_far_regret"],
            "final_top1_hits": summary[
                "final"
            ]["cumulative_top1_hits"],
            "first_top1_measurement": summary[
                "first_top1_measurement"
            ],
            "gate_opened_at_measurement": summary[
                "gate_opened_at_measurement"
            ],
            "assembly_active_rounds": summary[
                "assembly_active_rounds"
            ],
            "assembly_selected_count": summary[
                "assembly_selected_count"
            ],
            "assembly_outside_reservoir_count": summary[
                "assembly_outside_reservoir_count"
            ],
            "assembled_top1_hits": summary[
                "assembled_top1_hits"
            ],
            "selection_transcript_sha256": summary[
                "selection_transcript_sha256"
            ],
        }
        rows.append(row)
        by_landscape.setdefault(
            summary["landscape"],
            {},
        )[summary["condition"]] = summary

    comparison = pd.DataFrame(rows)
    comparison.to_csv(
        out / "PHASE2C_COMPARISON.csv",
        index=False,
    )

    verdict_rows = []
    wins = 0
    degradations = 0
    for landscape_name in (
        "GB1",
        "TRPB",
        "PHOQ",
    ):
        acquisition = by_landscape[
            landscape_name
        ]["acquisition_only"]
        combined = by_landscape[
            landscape_name
        ][
            "acquisition_then_gated_assembly"
        ]
        win = _landscape_win(
            acquisition,
            combined,
        )
        degradation = _joint_degradation(
            acquisition,
            combined,
        )
        wins += int(win)
        degradations += int(
            degradation
        )
        verdict_rows.append(
            {
                "landscape": landscape_name,
                "combined_win": bool(win),
                "joint_degradation": bool(
                    degradation
                ),
                "delta_auc": float(
                    combined[
                        "post_bootstrap_normalized_discovery_auc"
                    ]
                    - acquisition[
                        "post_bootstrap_normalized_discovery_auc"
                    ]
                ),
                "delta_regret": float(
                    combined["final"][
                        "best_so_far_regret"
                    ]
                    - acquisition["final"][
                        "best_so_far_regret"
                    ]
                ),
                "delta_top1_hits": int(
                    combined["final"][
                        "cumulative_top1_hits"
                    ]
                    - acquisition["final"][
                        "cumulative_top1_hits"
                    ]
                ),
            }
        )

    if wins >= 2 and degradations == 0:
        verdict = "PASS"
    elif wins >= 1:
        verdict = "PASS_WITH_LIMITATION"
    else:
        verdict = "NO_SCIENTIFIC_WIN"

    pd.DataFrame(
        verdict_rows
    ).to_csv(
        out / "PHASE2C_VERDICT_TABLE.csv",
        index=False,
    )

    manifest = {
        "version": "NABU_PHASE2C_MULTILANDSCAPE_V1",
        "scientific_claim": False,
        "stage": "2C_ASSEMBLY_PLUS_ACQUISITION_DEVELOPMENT_BENCHMARK",
        "verdict": verdict,
        "landscapes": landscape_records,
        "conditions": list(CONDITIONS),
        "reservoir_size": int(
            args.reservoir_size
        ),
        "bootstrap_count": int(
            args.bootstrap_count
        ),
        "target_count": int(
            args.target_count
        ),
        "batch_size": int(
            args.batch_size
        ),
        "assembly": {
            "target_order": 3,
            "beam_width": int(
                args.beam_width
            ),
            "proposal_frontier": int(
                args.proposal_frontier
            ),
            "min_support": int(
                args.min_support
            ),
        },
        "deterministic_replay": replay_checks,
        "verdict_table": verdict_rows,
        "comparison": rows,
        "blind_phase2d_data_loaded": False,
        "interpretation_rule": (
            "This is development evidence only. "
            "Final external winner claims require Phase 2D."
        ),
    }
    (
        out
        / "PHASE2C_MANIFEST.json"
    ).write_text(
        json.dumps(
            manifest,
            indent=2,
            default=_json_default,
        ),
        encoding="utf-8",
    )

    print(
        json.dumps(
            manifest,
            indent=2,
            default=_json_default,
        )
    )


if __name__ == "__main__":
    main()
