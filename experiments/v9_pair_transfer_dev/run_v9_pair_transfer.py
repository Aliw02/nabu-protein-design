from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import ndcg_score

from nabu_protein.higher_order import additive_score, model_diagnostics
from nabu_protein.v83 import NabuV83Model

SEED = 161
EXPECTED_IRED_SHA256 = (
    "aa45a2f85fb1af87b6b0e86397b3f8f292a061dc574be2536673b5bd490b0e74"
)
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def derive_reference(train_sequences: list[str]) -> str:
    if not train_sequences:
        raise ValueError("Empty fit/train identity set.")
    lengths = {len(seq) for seq in train_sequences}
    if len(lengths) != 1:
        raise ValueError(f"Fit/train sequences are not fixed-length: {lengths}")
    if any(not set(seq).issubset(STANDARD_AA) for seq in train_sequences):
        raise ValueError("Fit/train contains non-standard amino acids.")
    length = next(iter(lengths))
    return "".join(
        Counter(seq[i] for seq in train_sequences).most_common(1)[0][0]
        for i in range(length)
    )


def mutation_set(sequence: str, reference: str) -> tuple[str, ...]:
    sequence = str(sequence).strip().upper()
    if len(sequence) != len(reference):
        raise ValueError("Sequence/reference length mismatch.")
    if not set(sequence).issubset(STANDARD_AA):
        raise ValueError("Non-standard amino acid in sequence.")
    return tuple(
        f"{source}{i + 1}{target}"
        for i, (source, target) in enumerate(zip(reference, sequence))
        if source != target
    )


def evaluate(target: np.ndarray, prediction: np.ndarray) -> dict:
    target = np.asarray(target, dtype=float)
    prediction = np.asarray(prediction, dtype=float)
    if len(target) != len(prediction) or len(target) == 0:
        raise ValueError("Target/prediction length mismatch or empty evaluation set.")
    if not np.isfinite(target).all() or not np.isfinite(prediction).all():
        raise ValueError("Metrics require finite values.")

    rho = float(spearmanr(target, prediction).statistic)
    relevance = target - float(np.min(target))
    ndcg = float(ndcg_score(relevance[None, :], prediction[None, :]))

    n = len(target)
    k = int(math.ceil(0.01 * n))
    predicted_top = np.argsort(-prediction, kind="stable")[:k]
    true_top = np.argsort(-target, kind="stable")[:k]
    hits = len(set(predicted_top.tolist()).intersection(true_top.tolist()))

    best_true = float(np.max(target[predicted_top]))
    global_min = float(np.min(target))
    global_max = float(np.max(target))
    span = global_max - global_min
    regret = 0.0 if span <= 0.0 else float((global_max - best_true) / span)

    return {
        "spearman": rho,
        "ndcg": ndcg,
        "top1_percent_k": int(k),
        "top1_percent_hits": int(hits),
        "top1_percent_recall": float(hits / k),
        "top1_percent_enrichment": float(hits * n / (k * k)),
        "normalized_regret_top1pct": regret,
        "best_true_in_predicted_top1pct": best_true,
    }


def build_pair_factor_model(
    mutation_sets: list[tuple[str, ...]],
    labels: np.ndarray,
    oof_b2: np.ndarray,
    permute_targets: bool = False,
) -> dict:
    rows = [
        i
        for i, mutations in enumerate(mutation_sets)
        if len(mutations) == 2
    ]
    if not rows:
        raise RuntimeError("No double-mutant rows available for pair factorization.")

    pairs = [mutation_sets[i] for i in rows]
    residual = np.asarray(
        [float(labels[i] - oof_b2[i]) for i in rows],
        dtype=float,
    )

    if permute_targets:
        residual = np.random.default_rng(SEED).permutation(residual)

    nodes = sorted({mutation for pair in pairs for mutation in pair})
    node_index = {mutation: index for index, mutation in enumerate(nodes)}

    design = np.zeros((len(pairs), 1 + len(nodes)), dtype=float)
    design[:, 0] = 1.0
    support = Counter()

    for row_index, pair in enumerate(pairs):
        for mutation in pair:
            design[row_index, 1 + node_index[mutation]] = 1.0
            support[mutation] += 1

    coef, _, rank, singular_values = np.linalg.lstsq(
        design,
        residual,
        rcond=None,
    )

    return {
        "bias": float(coef[0]),
        "node_effect": {
            mutation: float(coef[1 + node_index[mutation]])
            for mutation in nodes
        },
        "node_support": {
            mutation: int(support[mutation])
            for mutation in nodes
        },
        "rows": int(len(pairs)),
        "nodes": int(len(nodes)),
        "rank": int(rank),
        "singular_min": (
            float(np.min(singular_values))
            if len(singular_values)
            else 0.0
        ),
        "target_mean": float(np.mean(residual)),
        "target_std": float(np.std(residual)),
        "permuted": bool(permute_targets),
    }


def transferred_pair_effect(pair: tuple[str, str], factor_model: dict) -> float | None:
    left, right = pair
    effect = factor_model["node_effect"]
    support = factor_model["node_support"]
    if left not in effect or right not in effect:
        return None

    left_conf = support[left] / (support[left] + 1.0)
    right_conf = support[right] / (support[right] + 1.0)
    confidence = math.sqrt(left_conf * right_conf)

    return float(
        (
            factor_model["bias"]
            + effect[left]
            + effect[right]
        )
        * confidence
    )


def score_v9_candidate(
    mutations: tuple[str, ...],
    base_model: dict,
    factor_model: dict,
) -> tuple[float, dict]:
    main_memory = base_model["main"]
    pair_memory = base_model["pair"]
    global_mean = float(base_model["global_mean"])

    if not all(mutation in main_memory for mutation in mutations):
        return global_mean, {
            "all_main_supported": False,
            "exact_pairs": 0,
            "transferred_pairs": 0,
            "unresolved_pairs": math.comb(len(mutations), 2)
            if len(mutations) >= 2
            else 0,
        }

    score = float(
        additive_score(
            mutations,
            global_mean,
            main_memory,
        )
    )
    exact_pairs = 0
    transferred_pairs = 0
    unresolved_pairs = 0

    for pair in combinations(mutations, 2):
        entry = pair_memory.get(pair)
        if entry is not None:
            score += float(entry["mean"] * entry["confidence"])
            exact_pairs += 1
            continue

        transfer = transferred_pair_effect(pair, factor_model)
        if transfer is None:
            unresolved_pairs += 1
        else:
            score += transfer
            transferred_pairs += 1

    return float(score), {
        "all_main_supported": True,
        "exact_pairs": int(exact_pairs),
        "transferred_pairs": int(transferred_pairs),
        "unresolved_pairs": int(unresolved_pairs),
    }


def b2_prediction(mutations: tuple[str, ...], base_model: dict) -> float:
    if not all(mutation in base_model["main"] for mutation in mutations):
        return float(base_model["global_mean"])
    return float(
        additive_score(
            mutations,
            float(base_model["global_mean"]),
            base_model["main"],
        )
    )


def phase2_abstention_prediction(
    model: NabuV83Model,
    mutation_sets: list[tuple[str, ...]],
    candidate_ids: list[str],
) -> tuple[np.ndarray, dict]:
    base = model.model["base"]
    strict_positions = []
    strict_mutations = []
    strict_ids = []

    for index, mutations in enumerate(mutation_sets):
        all_main = all(mutation in base["main"] for mutation in mutations)
        any_pair = any(
            pair in base["pair"]
            for pair in combinations(mutations, 2)
        )
        if all_main and any_pair:
            strict_positions.append(index)
            strict_mutations.append(mutations)
            strict_ids.append(candidate_ids[index])

    if not strict_ids:
        raise RuntimeError("No strict V8.3-scoreable test rows.")

    scored = model.score_candidates(strict_mutations, strict_ids)
    strict_scores = scored["V8_3_ADAPTIVE_ROUTER"].to_numpy(dtype=float)
    if not np.isfinite(strict_scores).all():
        raise RuntimeError("Frozen V8.3 router returned non-finite scores.")

    floor = float(np.min(strict_scores) - 1.0)
    prediction = np.full(len(mutation_sets), floor, dtype=float)
    for index, value in zip(strict_positions, strict_scores):
        prediction[index] = float(value)

    return prediction, {
        "strict_scoreable": int(len(strict_positions)),
        "abstained": int(len(mutation_sets) - len(strict_positions)),
        "abstention_floor": floor,
    }


def summarize_transfer_diagnostics(rows: list[dict]) -> dict:
    return {
        "all_main_supported_rows": int(
            sum(bool(row["all_main_supported"]) for row in rows)
        ),
        "rows_with_exact_pair": int(
            sum(int(row["exact_pairs"]) > 0 for row in rows)
        ),
        "rows_with_transferred_pair": int(
            sum(int(row["transferred_pairs"]) > 0 for row in rows)
        ),
        "exact_pair_contributions": int(
            sum(int(row["exact_pairs"]) for row in rows)
        ),
        "transferred_pair_contributions": int(
            sum(int(row["transferred_pairs"]) for row in rows)
        ),
        "unresolved_pair_contributions": int(
            sum(int(row["unresolved_pairs"]) for row in rows)
        ),
    }


def metrics_by_mutation_count(
    mutation_sets: list[tuple[str, ...]],
    target: np.ndarray,
    prediction: np.ndarray,
) -> dict:
    counts = np.asarray([len(value) for value in mutation_sets], dtype=int)
    result = {}
    for count in sorted(set(counts.tolist())):
        mask = counts == count
        if int(mask.sum()) < 10:
            continue
        result[str(count)] = evaluate(target[mask], prediction[mask])
    return result


def run(source_gz: Path, output_dir: Path) -> dict:
    actual_sha = sha256_file(source_gz)
    if actual_sha != EXPECTED_IRED_SHA256:
        raise RuntimeError(
            f"IRED SHA256 mismatch: expected {EXPECTED_IRED_SHA256}, got {actual_sha}"
        )

    frame = pd.read_csv(source_gz, compression="gzip")
    required = {"sequence", "target", "set", "validation"}
    missing = required - set(frame.columns)
    if missing:
        raise RuntimeError(f"IRED source missing columns: {sorted(missing)}")

    frame["sequence"] = frame["sequence"].astype(str).str.strip().str.upper()
    frame["target"] = pd.to_numeric(frame["target"], errors="raise")

    set_text = frame["set"].astype(str).str.lower()
    validation_true = frame["validation"].fillna(False).astype(bool)
    fit_frame = frame[set_text.eq("train") & ~validation_true].copy()
    validation_frame = frame[validation_true].copy()
    test_frame = frame[set_text.eq("test")].copy()

    if (len(fit_frame), len(validation_frame), len(test_frame)) != (3746, 662, 4178):
        raise RuntimeError(
            "Unexpected IRED split counts: "
            f"{len(fit_frame)}/{len(validation_frame)}/{len(test_frame)}"
        )

    reference = derive_reference(fit_frame["sequence"].tolist())
    fit_mutations = [
        mutation_set(sequence, reference)
        for sequence in fit_frame["sequence"]
    ]
    validation_mutations = [
        mutation_set(sequence, reference)
        for sequence in validation_frame["sequence"]
    ]
    test_mutations = [
        mutation_set(sequence, reference)
        for sequence in test_frame["sequence"]
    ]

    fit_labels = fit_frame["target"].to_numpy(dtype=float)
    validation_target = validation_frame["target"].to_numpy(dtype=float)
    test_target = test_frame["target"].to_numpy(dtype=float)

    model = NabuV83Model().fit(
        mutation_sets=fit_mutations,
        labels=fit_labels,
        candidate_ids=fit_frame["sequence"].tolist(),
    )
    base = model.model["base"]

    factor = build_pair_factor_model(
        fit_mutations,
        fit_labels,
        model.model["oof_b2"],
        permute_targets=False,
    )
    factor_perm = build_pair_factor_model(
        fit_mutations,
        fit_labels,
        model.model["oof_b2"],
        permute_targets=True,
    )

    phase2_prediction, phase2_diag = phase2_abstention_prediction(
        model,
        test_mutations,
        test_frame["sequence"].tolist(),
    )

    b2_test = np.asarray(
        [b2_prediction(mutations, base) for mutations in test_mutations],
        dtype=float,
    )

    v9_rows = [
        score_v9_candidate(mutations, base, factor)
        for mutations in test_mutations
    ]
    v9_test = np.asarray([row[0] for row in v9_rows], dtype=float)
    v9_diag_rows = [row[1] for row in v9_rows]

    perm_rows = [
        score_v9_candidate(mutations, base, factor_perm)
        for mutations in test_mutations
    ]
    perm_test = np.asarray([row[0] for row in perm_rows], dtype=float)

    validation_b2 = np.asarray(
        [b2_prediction(mutations, base) for mutations in validation_mutations],
        dtype=float,
    )
    validation_v9 = np.asarray(
        [
            score_v9_candidate(mutations, base, factor)[0]
            for mutations in validation_mutations
        ],
        dtype=float,
    )

    arms = {
        "PHASE2_V8_3_ABSTENTION": phase2_prediction,
        "B2_MAIN_ONLY": b2_test,
        "V9_PAIR_TRANSFER": v9_test,
        "V9_PAIR_TRANSFER_PERMUTED_CONTROL": perm_test,
    }
    metrics = {
        name: evaluate(test_target, prediction)
        for name, prediction in arms.items()
    }

    promising = bool(
        metrics["V9_PAIR_TRANSFER"]["spearman"]
        > metrics["B2_MAIN_ONLY"]["spearman"]
        and (
            metrics["V9_PAIR_TRANSFER"]["top1_percent_hits"]
            > metrics["B2_MAIN_ONLY"]["top1_percent_hits"]
            or metrics["V9_PAIR_TRANSFER"]["top1_percent_enrichment"]
            > metrics["B2_MAIN_ONLY"]["top1_percent_enrichment"]
        )
        and metrics["V9_PAIR_TRANSFER"]["normalized_regret_top1pct"]
        <= metrics["B2_MAIN_ONLY"]["normalized_regret_top1pct"]
        and metrics["V9_PAIR_TRANSFER"]["spearman"]
        > metrics["V9_PAIR_TRANSFER_PERMUTED_CONTROL"]["spearman"]
        and np.isfinite(v9_test).all()
    )

    predictions = pd.DataFrame(
        {
            "candidate_id": test_frame["sequence"].tolist(),
            "mutation_count": [len(value) for value in test_mutations],
            "target": test_target,
            **arms,
            "all_main_supported": [
                row["all_main_supported"] for row in v9_diag_rows
            ],
            "exact_pairs": [row["exact_pairs"] for row in v9_diag_rows],
            "transferred_pairs": [
                row["transferred_pairs"] for row in v9_diag_rows
            ],
            "unresolved_pairs": [
                row["unresolved_pairs"] for row in v9_diag_rows
            ],
        }
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    predictions.to_csv(output_dir / "V9_IRED_PREDICTIONS.csv", index=False)

    result = {
        "version": "NABU_V9_0_PAIR_TRANSFER_DEV_V1",
        "status": (
            "MECHANISTICALLY_PROMISING"
            if promising
            else "DEVELOPMENT_FAIL_OR_INCONCLUSIVE"
        ),
        "scientific_claim": False,
        "phase3_open": False,
        "source_sha256": actual_sha,
        "base_frozen_v83_commit": (
            "c8afdcd7698231d95a39ef13a3fe22b6e3f507c5"
        ),
        "split_counts": {
            "fit": int(len(fit_frame)),
            "validation_held_out": int(len(validation_frame)),
            "test_post_blind_development": int(len(test_frame)),
        },
        "reference_sha256": hashlib.sha256(reference.encode("utf-8")).hexdigest(),
        "reference_length": int(len(reference)),
        "router_decision": model.router_decision,
        "v83_model_diagnostics": model_diagnostics(model.model),
        "phase2_abstention_reproduction": phase2_diag,
        "factor_model": {
            key: value
            for key, value in factor.items()
            if key not in {"node_effect", "node_support"}
        },
        "factor_node_count": int(len(factor["node_effect"])),
        "factor_node_support_min": int(min(factor["node_support"].values())),
        "factor_node_support_max": int(max(factor["node_support"].values())),
        "transfer_diagnostics": summarize_transfer_diagnostics(v9_diag_rows),
        "validation_diagnostic": {
            "B2_MAIN_ONLY": evaluate(validation_target, validation_b2),
            "V9_PAIR_TRANSFER": evaluate(validation_target, validation_v9),
            "role": "HELD_OUT_DIAGNOSTIC_NOT_USED_FOR_SELECTION",
        },
        "test_metrics": metrics,
        "deltas": {
            "V9_minus_B2_spearman": float(
                metrics["V9_PAIR_TRANSFER"]["spearman"]
                - metrics["B2_MAIN_ONLY"]["spearman"]
            ),
            "V9_minus_permuted_spearman": float(
                metrics["V9_PAIR_TRANSFER"]["spearman"]
                - metrics["V9_PAIR_TRANSFER_PERMUTED_CONTROL"]["spearman"]
            ),
            "V9_minus_B2_top1_hits": int(
                metrics["V9_PAIR_TRANSFER"]["top1_percent_hits"]
                - metrics["B2_MAIN_ONLY"]["top1_percent_hits"]
            ),
            "V9_minus_B2_regret": float(
                metrics["V9_PAIR_TRANSFER"]["normalized_regret_top1pct"]
                - metrics["B2_MAIN_ONLY"]["normalized_regret_top1pct"]
            ),
        },
        "v9_metrics_by_mutation_count": metrics_by_mutation_count(
            test_mutations,
            test_target,
            v9_test,
        ),
        "decision_rule_pass": promising,
        "interpretation_boundary": (
            "IRED targets were already revealed by Phase 2D. "
            "This is post-blind V9 development evidence only and "
            "cannot change the Phase-2 verdict."
        ),
    }

    (output_dir / "V9_IRED_RESULTS.json").write_text(
        json.dumps(result, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2))
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source_gz")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    run(Path(args.source_gz), Path(args.out))


if __name__ == "__main__":
    main()
