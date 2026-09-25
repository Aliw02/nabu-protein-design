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
from nabu_protein.v83 import NabuV83Model  # noqa: E402


def sha256_ids(candidate_ids: list[str]) -> str:
    payload = "\n".join(candidate_ids).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def extract_vocabulary(identity_pool: pd.DataFrame) -> list[str]:
    tokens = sorted(
        {
            token
            for mutation_set in identity_pool["mutation_set"]
            for token in mutation_set
        },
        key=lambda token: (int(token[1:-1]), token),
    )
    return tokens


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run Phase-2B assembly-only RhlA development POC."
    )
    parser.add_argument(
        "--input",
        default="rhla_sample_efficiency_sealed_input/TRAINING_POOL.csv",
    )
    parser.add_argument(
        "--out",
        default="phase2/2B_candidate_assembly/results/rhla_assembly_v1",
    )
    parser.add_argument("--seed-fraction", type=float, default=0.05)
    parser.add_argument("--target-order", type=int, default=3)
    parser.add_argument("--beam-width", type=int, default=256)
    parser.add_argument("--proposal-count", type=int, default=64)
    args = parser.parse_args()

    source = pd.read_csv(args.input)
    identity_pool = identity_pool_from_frame(source)
    truth = truth_frame_from_frame(source)

    vocabulary_tokens = extract_vocabulary(identity_pool)
    reference_sequence = infer_reference_sequence_from_tokens(
        vocabulary_tokens
    )
    reference = ReferenceProtein(
        name="RhlA_development_reference_inferred_from_vocabulary",
        sequence=reference_sequence,
    )
    vocabulary = MutationVocabulary.from_tokens(
        reference,
        vocabulary_tokens,
    )

    seed_count = int(math.ceil(len(identity_pool) * args.seed_fraction))
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

    if len(measured) != seed_count:
        raise RuntimeError("Seed reveal count mismatch.")

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
        target_order=args.target_order,
        beam_width=args.beam_width,
        proposal_count=args.proposal_count,
    )

    proposal_ids = proposals["candidate_id"].astype(str).tolist()
    proposal_hash = sha256_ids(proposal_ids)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    pre_eval = {
        "stage": "PROPOSALS_FROZEN_BEFORE_VIRTUAL_TRUTH_LOOKUP",
        "proposal_count": int(len(proposal_ids)),
        "proposal_ids": proposal_ids,
        "proposal_sha256": proposal_hash,
        "target_order": int(args.target_order),
        "beam_width": int(args.beam_width),
        "reference_name": reference.name,
        "reference_sequence": reference.sequence,
        "vocabulary_size": int(len(vocabulary.tokens)),
        "truth_exposed_to_assembler_before_freeze": False,
        "assembly_trace": trace,
    }
    (out / "PROPOSAL_FREEZE.json").write_text(
        json.dumps(pre_eval, indent=2),
        encoding="utf-8",
    )

    proposal_frame = proposals.copy()
    proposal_frame["mutation_set"] = proposal_frame[
        "mutation_set"
    ].map(lambda value: ":".join(value))
    proposal_frame.to_csv(out / "PROPOSALS.csv", index=False)

    source_truth = source[["candidate_id", "DMS_score"]].copy()
    source_truth["candidate_id"] = source_truth["candidate_id"].astype(str)
    source_truth["DMS_score"] = source_truth["DMS_score"].astype(float)

    evaluated = proposals[["candidate_id"]].merge(
        source_truth,
        on="candidate_id",
        how="left",
        validate="one_to_one",
    )
    matched = evaluated["DMS_score"].notna()
    matched_values = evaluated.loc[matched, "DMS_score"].to_numpy(dtype=float)

    full_values = source_truth["DMS_score"].to_numpy(dtype=float)
    top1_cutoff = float(np.quantile(full_values, 0.99))

    summary = {
        "version": "NABU_PHASE2B_ASSEMBLY_POC_V1",
        "scientific_claim": False,
        "stage": "2B_ASSEMBLY_ONLY_DEVELOPMENT",
        "input": str(args.input),
        "source_pool_size": int(len(source)),
        "bootstrap_count": int(seed_count),
        "reference_sequence": reference.sequence,
        "vocabulary_size": int(len(vocabulary.tokens)),
        "target_order": int(args.target_order),
        "beam_width": int(args.beam_width),
        "proposal_count": int(len(proposal_ids)),
        "proposal_sha256_before_truth_lookup": proposal_hash,
        "proposal_scoreable_count": int(proposals["scoreable"].sum()),
        "proposal_measured_overlap": int(
            len(set(proposal_ids) & seed_set)
        ),
        "proposal_duplicate_count": int(
            len(proposal_ids) - len(set(proposal_ids))
        ),
        "oracle_match_count": int(matched.sum()),
        "oracle_unmatched_count": int((~matched).sum()),
        "matched_best_true_fitness": (
            None
            if len(matched_values) == 0
            else float(np.max(matched_values))
        ),
        "matched_top1_hits": int(
            np.sum(matched_values >= top1_cutoff)
        ),
        "source_top1_cutoff": top1_cutoff,
        "assembly_trace": trace,
        "interpretation_rule": (
            "Historical RhlA truth is consulted only after proposal freeze. "
            "Unmatched proposals are novel-to-table proposals with unknown "
            "virtual truth and are not failures."
        ),
    }

    (out / "ASSEMBLY_SUMMARY.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
