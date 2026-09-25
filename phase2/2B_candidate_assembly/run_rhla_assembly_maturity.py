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
PHASE2A0 = HERE.parent / "2A_closed_loop_acquisition" / "2A0_campaign_simulator"
sys.path.insert(0, str(PHASE2A0))

from campaign import (  # noqa: E402
    LABEL_COLUMN,
    bootstrap_seed_ids,
    identity_pool_from_frame,
    truth_frame_from_frame,
)
from assembly import (  # noqa: E402
    CandidateAssembler,
    MutationVocabulary,
    ReferenceProtein,
    infer_reference_sequence_from_tokens,
)
from nabu_protein.higher_order import model_diagnostics  # noqa: E402
from nabu_protein.v83 import NabuV83Model  # noqa: E402


def sha256_ids(candidate_ids: list[str]) -> str:
    return hashlib.sha256(
        "\n".join(candidate_ids).encode("utf-8")
    ).hexdigest()


def extract_vocabulary(identity_pool: pd.DataFrame) -> list[str]:
    return sorted(
        {
            token
            for mutation_set in identity_pool["mutation_set"]
            for token in mutation_set
        },
        key=lambda token: (int(token[1:-1]), token),
    )


def percentile(values: np.ndarray, value: float) -> float:
    values = np.asarray(values, dtype=float)
    return float(
        (
            np.sum(values < value)
            + 0.5 * np.sum(values == value)
        )
        / len(values)
    )


def run_level(
    fraction: float,
    source: pd.DataFrame,
    identity_pool: pd.DataFrame,
    truth: pd.DataFrame,
    vocabulary: MutationVocabulary,
    out: Path,
    target_order: int,
    beam_width: int,
    proposal_count: int,
) -> dict:
    seed_count = int(math.ceil(len(identity_pool) * fraction))
    seed_ids = bootstrap_seed_ids(
        identity_pool,
        target_count=seed_count,
        min_per_fold=1,
    )
    seed_set = set(seed_ids)

    seed_identity = identity_pool[
        identity_pool["candidate_id"].isin(seed_set)
    ].copy()
    seed_truth = truth[
        truth["candidate_id"].isin(seed_set)
    ].copy()
    measured = seed_identity.merge(
        seed_truth,
        on="candidate_id",
        how="inner",
        validate="one_to_one",
    )

    model = NabuV83Model().fit(
        measured["mutation_set"].tolist(),
        measured[LABEL_COLUMN].to_numpy(dtype=float),
        measured["candidate_id"].astype(str).tolist(),
    )

    assembler = CandidateAssembler(
        model=model,
        vocabulary=vocabulary,
        measured_mutation_sets=measured["mutation_set"].tolist(),
    )
    proposals, trace = assembler.assemble(
        target_order=target_order,
        beam_width=beam_width,
        proposal_count=proposal_count,
    )

    proposal_ids = proposals["candidate_id"].astype(str).tolist()
    proposal_hash = sha256_ids(proposal_ids)

    source_truth = source[["candidate_id", "DMS_score"]].copy()
    source_truth["candidate_id"] = source_truth["candidate_id"].astype(str)
    source_truth["DMS_score"] = source_truth["DMS_score"].astype(float)

    evaluated = proposals[["candidate_id"]].merge(
        source_truth,
        on="candidate_id",
        how="left",
        validate="one_to_one",
    )
    matched_values = evaluated.loc[
        evaluated["DMS_score"].notna(),
        "DMS_score",
    ].to_numpy(dtype=float)

    full_values = source_truth["DMS_score"].to_numpy(dtype=float)
    top1_cutoff = float(np.quantile(full_values, 0.99))
    matched_percentiles = np.array(
        [percentile(full_values, value) for value in matched_values],
        dtype=float,
    )

    diagnostics = model_diagnostics(model.model)

    level_dir = out / f"evidence_{int(round(fraction * 100)):02d}pct"
    level_dir.mkdir(parents=True, exist_ok=True)

    proposal_frame = proposals.copy()
    proposal_frame["mutation_set"] = proposal_frame[
        "mutation_set"
    ].map(lambda value: ":".join(value))
    proposal_frame.to_csv(
        level_dir / "PROPOSALS.csv",
        index=False,
    )

    summary = {
        "fraction": float(fraction),
        "seed_count": int(seed_count),
        "proposal_count": int(len(proposal_ids)),
        "proposal_sha256_before_truth_lookup": proposal_hash,
        "proposal_scoreable_count": int(proposals["scoreable"].sum()),
        "proposal_measured_overlap": int(
            len(set(proposal_ids) & seed_set)
        ),
        "oracle_match_count": int(len(matched_values)),
        "oracle_unmatched_count": int(
            len(proposal_ids) - len(matched_values)
        ),
        "matched_best_true_fitness": (
            None
            if len(matched_values) == 0
            else float(np.max(matched_values))
        ),
        "matched_top1_hits": int(
            np.sum(matched_values >= top1_cutoff)
        ),
        "matched_mean_percentile": (
            None
            if len(matched_percentiles) == 0
            else float(np.mean(matched_percentiles))
        ),
        "matched_median_percentile": (
            None
            if len(matched_percentiles) == 0
            else float(np.median(matched_percentiles))
        ),
        "router_mode": str(model.router_decision["mode"]),
        "model_diagnostics": diagnostics,
        "assembly_trace": trace,
    }
    (level_dir / "SUMMARY.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run RhlA assembly evidence-maturity ablation."
    )
    parser.add_argument(
        "--input",
        default="rhla_sample_efficiency_sealed_input/TRAINING_POOL.csv",
    )
    parser.add_argument(
        "--out",
        default=(
            "phase2/2B_candidate_assembly/results/"
            "rhla_assembly_maturity_v1"
        ),
    )
    parser.add_argument("--target-order", type=int, default=3)
    parser.add_argument("--beam-width", type=int, default=256)
    parser.add_argument("--proposal-count", type=int, default=64)
    args = parser.parse_args()

    source = pd.read_csv(args.input)
    identity_pool = identity_pool_from_frame(source)
    truth = truth_frame_from_frame(source)

    vocabulary_tokens = extract_vocabulary(identity_pool)
    reference = ReferenceProtein(
        name="RhlA_development_reference_inferred_from_vocabulary",
        sequence=infer_reference_sequence_from_tokens(vocabulary_tokens),
    )
    vocabulary = MutationVocabulary.from_tokens(
        reference,
        vocabulary_tokens,
    )

    fractions = [0.05, 0.10, 0.20]

    seed_sets = []
    for fraction in fractions:
        count = int(math.ceil(len(identity_pool) * fraction))
        seed_sets.append(
            set(
                bootstrap_seed_ids(
                    identity_pool,
                    target_count=count,
                    min_per_fold=1,
                )
            )
        )
    if not (
        seed_sets[0].issubset(seed_sets[1])
        and seed_sets[1].issubset(seed_sets[2])
    ):
        raise RuntimeError("Evidence-level seed memberships are not nested.")

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    summaries = []
    for fraction in fractions:
        print(f"[2B-MATURITY] START {fraction:.2f}")
        summary = run_level(
            fraction=fraction,
            source=source,
            identity_pool=identity_pool,
            truth=truth,
            vocabulary=vocabulary,
            out=out,
            target_order=args.target_order,
            beam_width=args.beam_width,
            proposal_count=args.proposal_count,
        )
        summaries.append(summary)
        print(
            f"[2B-MATURITY] DONE {fraction:.2f} "
            f"matched={summary['oracle_match_count']} "
            f"best={summary['matched_best_true_fitness']} "
            f"top1={summary['matched_top1_hits']} "
            f"mean_pct={summary['matched_mean_percentile']}"
        )

    rows = []
    for summary in summaries:
        diagnostics = summary["model_diagnostics"]
        rows.append(
            {
                "fraction": summary["fraction"],
                "seed_count": summary["seed_count"],
                "proposal_scoreable_count": summary[
                    "proposal_scoreable_count"
                ],
                "oracle_match_count": summary["oracle_match_count"],
                "matched_best_true_fitness": summary[
                    "matched_best_true_fitness"
                ],
                "matched_top1_hits": summary["matched_top1_hits"],
                "matched_mean_percentile": summary[
                    "matched_mean_percentile"
                ],
                "matched_median_percentile": summary[
                    "matched_median_percentile"
                ],
                "router_mode": summary["router_mode"],
                "pair_entries": diagnostics["pair"]["entries"],
                "triplet_entries": diagnostics[
                    "triplet_crossfit"
                ]["entries"],
                "quartet_entries": diagnostics[
                    "quartet_crossfit"
                ]["entries"],
            }
        )

    pd.DataFrame(rows).to_csv(
        out / "ASSEMBLY_MATURITY_COMPARISON.csv",
        index=False,
    )

    manifest = {
        "version": "NABU_PHASE2B_ASSEMBLY_MATURITY_V1",
        "scientific_claim": False,
        "stage": "2B_ASSEMBLY_EVIDENCE_MATURITY",
        "single_changed_variable": "measured_evidence_fraction",
        "fractions": fractions,
        "nested_identity_only_seeds": True,
        "target_order": int(args.target_order),
        "beam_width": int(args.beam_width),
        "proposal_count": int(args.proposal_count),
        "comparison": rows,
    }
    (out / "ASSEMBLY_MATURITY_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )

    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
