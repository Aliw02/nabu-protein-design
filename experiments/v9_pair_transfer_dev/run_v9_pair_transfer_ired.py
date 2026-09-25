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
EXPECTED_IRED_SHA256 = "aa45a2f85fb1af87b6b0e86397b3f8f292a061dc574be2536673b5bd490b0e74"
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
VERSION = "NABU_V9_0_PAIR_TRANSFER_DEV_V1"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def derive_reference(sequences: list[str]) -> str:
    if not sequences:
        raise ValueError("Empty fit sequence set.")
    lengths = {len(sequence) for sequence in sequences}
    if len(lengths) != 1:
        raise ValueError("Fit sequences are not fixed length.")
    if any(not set(sequence).issubset(STANDARD_AA) for sequence in sequences):
        raise ValueError("Fit sequences contain non-standard amino acids.")
    length = next(iter(lengths))
    return "".join(
        Counter(sequence[i] for sequence in sequences).most_common(1)[0][0]
        for i in range(length)
    )


def mutation_set(sequence: str, reference: str) -> tuple[str, ...]:
    sequence = str(sequence).strip().upper()
    if len(sequence) != len(reference):
        raise ValueError("Sequence/reference length mismatch.")
    if not set(sequence).issubset(STANDARD_AA):
        raise ValueError("Sequence contains non-standard amino acids.")
    return tuple(
        f"{source}{i + 1}{target}"
        for i, (source, target) in enumerate(zip(reference, sequence))
        if source != target
    )


def finite_spearman(target, prediction):
    value = float(spearmanr(target, prediction).statistic)
    return value if np.isfinite(value) else None


def evaluate(target, prediction) -> dict:
    target = np.asarray(target, dtype=float)
    prediction = np.asarray(prediction, dtype=float)
    if len(target) != len(prediction) or len(target) == 0:
        raise ValueError("Target/prediction mismatch.")
    if not np.isfinite(target).all() or not np.isfinite(prediction).all():
        raise ValueError("Evaluation requires finite values.")

    relevance = target - float(np.min(target))
    n = len(target)
    k = int(math.ceil(0.01 * n))
    predicted_top = np.argsort(-prediction, kind="stable")[:k]
    true_top = np.argsort(-target, kind="stable")[:k]
    hits = int(len(set(predicted_top.tolist()) & set(true_top.tolist())))
    best_true = float(np.max(target[predicted_top]))
    low = float(np.min(target))
    high = float(np.max(target))
    span = high - low

    return {
        "count": int(n),
        "spearman": finite_spearman(target, prediction),
        "ndcg": float(ndcg_score(relevance[None, :], prediction[None, :])),
        "top1_percent_k": int(k),
        "top1_percent_hits": hits,
        "top1_percent_recall": float(hits / k),
        "top1_percent_enrichment": float(hits * n / (k * k)),
        "normalized_regret_top1pct": 0.0 if span <= 0 else float((high - best_true) / span),
        "best_true_in_predicted_top1pct": best_true,
    }


def prediction_hash(candidate_ids, values) -> str:
    payload = "\n".join(
        f"{candidate_id},{float(value):.17g}"
        for candidate_id, value in zip(candidate_ids, values)
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def fit_pair_factorization(
    model,
    mutation_sets,
    labels,
    *,
    permute_targets: bool = False,
) -> dict:
    labels = np.asarray(labels, dtype=float)
    oof_b2 = np.asarray(model.model["oof_b2"], dtype=float)
    double_indices = np.array(
        [i for i, mutations in enumerate(mutation_sets) if len(mutations) == 2],
        dtype=int,
    )
    if len(double_indices) == 0:
        raise RuntimeError("No double-mutant rows available.")

    edges = [mutation_sets[i] for i in double_indices]
    residual = labels[double_indices] - oof_b2[double_indices]
    if permute_targets:
        residual = np.random.default_rng(SEED).permutation(residual)

    nodes = sorted({mutation for edge in edges for mutation in edge})
    node_index = {mutation: i for i, mutation in enumerate(nodes)}
    design = np.zeros((len(edges), len(nodes) + 1), dtype=float)
    design[:, 0] = 1.0
    support = Counter()

    for row, edge in enumerate(edges):
        for mutation in edge:
            design[row, 1 + node_index[mutation]] = 1.0
            support[mutation] += 1

    coefficients, _, rank, singular_values = np.linalg.lstsq(
        design,
        residual,
        rcond=None,
    )
    fitted = design @ coefficients
    error = residual - fitted

    return {
        "bias": float(coefficients[0]),
        "propensity": {
            mutation: float(coefficients[1 + node_index[mutation]])
            for mutation in nodes
        },
        "confidence": {
            mutation: float(support[mutation] / (support[mutation] + 1.0))
            for mutation in nodes
        },
        "support": {mutation: int(support[mutation]) for mutation in nodes},
        "double_count": int(len(edges)),
        "node_count": int(len(nodes)),
        "design_rank": int(rank),
        "singular_value_min": float(np.min(singular_values)),
        "singular_value_max": float(np.max(singular_values)),
        "fit_rmse": float(np.sqrt(np.mean(error ** 2))),
        "target_mean": float(np.mean(residual)),
        "target_std": float(np.std(residual)),
        "permuted": bool(permute_targets),
    }


def transferred_pair_effect(pair, factor):
    left, right = pair
    propensity = factor["propensity"]
    confidence = factor["confidence"]
    if left not in propensity or right not in propensity:
        return None
    raw = float(factor["bias"] + propensity[left] + propensity[right])
    weight = math.sqrt(float(confidence[left]) * float(confidence[right]))
    return float(raw * weight)


def score_v9(mutations, model, factor):
    base = model.model["base"]
    main_memory = base["main"]
    pair_memory = base["pair"]
    global_mean = float(base["global_mean"])

    if not all(mutation in main_memory for mutation in mutations):
        return global_mean, {
            "all_main_supported": False,
            "exact_pair_count": 0,
            "transferred_pair_count": 0,
            "unsupported_pair_count": int(math.comb(len(mutations), 2))
            if len(mutations) >= 2 else 0,
        }

    score = float(additive_score(mutations, global_mean, main_memory))
    exact_count = 0
    transferred_count = 0
    unsupported_count = 0

    for pair in combinations(mutations, 2):
        exact = pair_memory.get(pair)
        if exact is not None:
            score += float(exact["mean"] * exact["confidence"])
            exact_count += 1
            continue

        transferred = transferred_pair_effect(pair, factor)
        if transferred is None:
            unsupported_count += 1
        else:
            score += transferred
            transferred_count += 1

    return score, {
        "all_main_supported": True,
        "exact_pair_count": int(exact_count),
        "transferred_pair_count": int(transferred_count),
        "unsupported_pair_count": int(unsupported_count),
    }


def phase2_abstention_predictions(model, mutation_sets, candidate_ids):
    base = model.model["base"]
    main_memory = base["main"]
    pair_memory = base["pair"]
    scoreable = np.array(
        [
            all(mutation in main_memory for mutation in mutations)
            and any(pair in pair_memory for pair in combinations(mutations, 2))
            for mutations in mutation_sets
        ],
        dtype=bool,
    )
    prediction = np.full(len(mutation_sets), np.nan, dtype=float)
    if scoreable.any():
        indices = np.where(scoreable)[0]
        scored = model.score_candidates(
            [mutation_sets[i] for i in indices],
            [candidate_ids[i] for i in indices],
        )
        values = scored["V8_3_ADAPTIVE_ROUTER"].to_numpy(dtype=float)
        prediction[indices] = values
        floor = float(np.min(values) - 1.0)
    else:
        floor = float(base["global_mean"] - 1.0)
    prediction[~scoreable] = floor
    return prediction, scoreable


def b2_predictions(model, mutation_sets):
    base = model.model["base"]
    main_memory = base["main"]
    global_mean = float(base["global_mean"])
    all_main = np.array(
        [
            all(mutation in main_memory for mutation in mutations)
            for mutations in mutation_sets
        ],
        dtype=bool,
    )
    prediction = np.array(
        [
            float(additive_score(mutations, global_mean, main_memory))
            if supported else global_mean
            for mutations, supported in zip(mutation_sets, all_main)
        ],
        dtype=float,
    )
    return prediction, all_main


def load_ired(path: Path):
    actual_sha = sha256_file(path)
    if actual_sha != EXPECTED_IRED_SHA256:
        raise RuntimeError(
            f"IRED SHA256 mismatch: expected {EXPECTED_IRED_SHA256}, got {actual_sha}"
        )

    frame = pd.read_csv(path, compression="gzip")
    required = {"sequence", "target", "set", "validation"}
    missing = required - set(frame.columns)
    if missing:
        raise RuntimeError(f"IRED source missing columns: {sorted(missing)}")

    frame["sequence"] = frame["sequence"].astype(str).str.strip().str.upper()
    frame["target"] = pd.to_numeric(frame["target"], errors="raise")
    set_text = frame["set"].astype(str).str.lower()
    validation_true = frame["validation"].fillna(False).astype(bool)

    fit = frame[set_text.eq("train") & ~validation_true].copy()
    validation = frame[validation_true].copy()
    test = frame[set_text.eq("test")].copy()

    if (len(fit), len(validation), len(test)) != (3746, 662, 4178):
        raise RuntimeError(
            f"Unexpected IRED split counts: {len(fit)}/{len(validation)}/{len(test)}"
        )
    return actual_sha, fit, validation, test


def evaluate_by_mutation_count(target, prediction, mutation_sets):
    counts = np.array([len(mutations) for mutations in mutation_sets], dtype=int)
    result = {}
    for count in sorted(set(counts.tolist())):
        mask = counts == count
        if int(mask.sum()) >= 2:
            result[str(count)] = evaluate(target[mask], prediction[mask])
    return result


def run(source_gz: Path, output_dir: Path) -> dict:
    source_sha, fit_frame, validation_frame, test_frame = load_ired(source_gz)
    reference = derive_reference(fit_frame["sequence"].tolist())

    fit_mutations = [mutation_set(x, reference) for x in fit_frame["sequence"]]
    validation_mutations = [
        mutation_set(x, reference) for x in validation_frame["sequence"]
    ]
    test_mutations = [mutation_set(x, reference) for x in test_frame["sequence"]]

    fit_target = fit_frame["target"].to_numpy(dtype=float)
    validation_target = validation_frame["target"].to_numpy(dtype=float)
    test_target = test_frame["target"].to_numpy(dtype=float)

    model = NabuV83Model().fit(
        mutation_sets=fit_mutations,
        labels=fit_target,
        candidate_ids=fit_frame["sequence"].tolist(),
    )

    factor = fit_pair_factorization(
        model,
        fit_mutations,
        fit_target,
        permute_targets=False,
    )
    permuted_factor = fit_pair_factorization(
        model,
        fit_mutations,
        fit_target,
        permute_targets=True,
    )

    test_ids = test_frame["sequence"].tolist()
    phase2_prediction, strict_scoreable = phase2_abstention_predictions(
        model,
        test_mutations,
        test_ids,
    )
    b2_prediction, all_main = b2_predictions(model, test_mutations)

    v9_rows = [score_v9(mutations, model, factor) for mutations in test_mutations]
    v9_prediction = np.array([value for value, _ in v9_rows], dtype=float)
    v9_diag = [diag for _, diag in v9_rows]

    permuted_prediction = np.array(
        [
            score_v9(mutations, model, permuted_factor)[0]
            for mutations in test_mutations
        ],
        dtype=float,
    )

    arms = {
        "PHASE2_V8_3_ABSTENTION": phase2_prediction,
        "B2_MAIN_ONLY": b2_prediction,
        "V9_PAIR_TRANSFER": v9_prediction,
        "V9_PAIR_TRANSFER_PERMUTED_CONTROL": permuted_prediction,
    }
    if not all(np.isfinite(values).all() for values in arms.values()):
        raise RuntimeError("A development arm produced non-finite predictions.")

    metrics = {
        name: evaluate(test_target, values)
        for name, values in arms.items()
    }
    by_mutation_count = {
        name: evaluate_by_mutation_count(test_target, values, test_mutations)
        for name, values in arms.items()
    }

    non_strict_mask = all_main & ~strict_scoreable
    non_strict_metrics = {
        "count": int(non_strict_mask.sum()),
        "B2_MAIN_ONLY": evaluate(
            test_target[non_strict_mask],
            b2_prediction[non_strict_mask],
        ),
        "V9_PAIR_TRANSFER": evaluate(
            test_target[non_strict_mask],
            v9_prediction[non_strict_mask],
        ),
        "V9_PAIR_TRANSFER_PERMUTED_CONTROL": evaluate(
            test_target[non_strict_mask],
            permuted_prediction[non_strict_mask],
        ),
    }

    validation_b2, validation_all_main = b2_predictions(
        model,
        validation_mutations,
    )
    validation_v9 = np.array(
        [score_v9(mutations, model, factor)[0] for mutations in validation_mutations],
        dtype=float,
    )
    validation_metrics = {
        "all_main_supported": int(validation_all_main.sum()),
        "B2_MAIN_ONLY": evaluate(validation_target, validation_b2),
        "V9_PAIR_TRANSFER": evaluate(validation_target, validation_v9),
    }

    v9_spearman = metrics["V9_PAIR_TRANSFER"]["spearman"]
    b2_spearman = metrics["B2_MAIN_ONLY"]["spearman"]
    perm_spearman = metrics["V9_PAIR_TRANSFER_PERMUTED_CONTROL"]["spearman"]
    non_strict_v9 = non_strict_metrics["V9_PAIR_TRANSFER"]["spearman"]
    non_strict_b2 = non_strict_metrics["B2_MAIN_ONLY"]["spearman"]

    checks = {
        "spearman_beats_b2": bool(
            v9_spearman is not None
            and b2_spearman is not None
            and v9_spearman > b2_spearman
        ),
        "elite_metric_improves": bool(
            metrics["V9_PAIR_TRANSFER"]["top1_percent_hits"]
            > metrics["B2_MAIN_ONLY"]["top1_percent_hits"]
            or metrics["V9_PAIR_TRANSFER"]["top1_percent_enrichment"]
            > metrics["B2_MAIN_ONLY"]["top1_percent_enrichment"]
        ),
        "normalized_regret_not_worse": bool(
            metrics["V9_PAIR_TRANSFER"]["normalized_regret_top1pct"]
            <= metrics["B2_MAIN_ONLY"]["normalized_regret_top1pct"]
        ),
        "spearman_beats_permuted_control": bool(
            v9_spearman is not None
            and perm_spearman is not None
            and v9_spearman > perm_spearman
        ),
        "non_strict_region_beats_b2": bool(
            non_strict_v9 is not None
            and non_strict_b2 is not None
            and non_strict_v9 > non_strict_b2
        ),
        "finite_predictions": True,
    }

    diagnostics = pd.DataFrame(v9_diag)
    result = {
        "version": VERSION,
        "status": "DEVELOPMENT_DIAGNOSTIC_AFTER_PHASE2_REVEAL",
        "phase3_opened": False,
        "scientific_claim": False,
        "source_sha256": source_sha,
        "frozen_v83_core": "c8afdcd7698231d95a39ef13a3fe22b6e3f507c5",
        "split_counts": {
            "fit": int(len(fit_frame)),
            "validation": int(len(validation_frame)),
            "test": int(len(test_frame)),
        },
        "reference_sha256": hashlib.sha256(reference.encode("utf-8")).hexdigest(),
        "reference_length": int(len(reference)),
        "v83_router_decision": model.router_decision,
        "v83_model_diagnostics": model_diagnostics(model.model),
        "factorization": {
            key: value
            for key, value in factor.items()
            if key not in {"propensity", "confidence", "support"}
        },
        "permuted_factorization": {
            key: value
            for key, value in permuted_factor.items()
            if key not in {"propensity", "confidence", "support"}
        },
        "coverage": {
            "strict_v83_scoreable": int(strict_scoreable.sum()),
            "all_main_supported": int(all_main.sum()),
            "v9_non_strict_transfer_region": int(non_strict_mask.sum()),
            "exact_pair_contributions_total": int(
                diagnostics["exact_pair_count"].sum()
            ),
            "transferred_pair_contributions_total": int(
                diagnostics["transferred_pair_count"].sum()
            ),
            "unsupported_pair_contributions_total": int(
                diagnostics["unsupported_pair_count"].sum()
            ),
        },
        "metrics": metrics,
        "non_strict_transfer_region_metrics": non_strict_metrics,
        "validation_diagnostic": validation_metrics,
        "by_mutation_count": by_mutation_count,
        "prediction_sha256": {
            name: prediction_hash(test_ids, values)
            for name, values in arms.items()
        },
        "decision_checks": checks,
        "development_decision": (
            "V9_PAIR_TRANSFER_PROMISING"
            if all(checks.values())
            else "V9_PAIR_TRANSFER_NOT_YET_SUPPORTED"
        ),
        "interpretation_boundary": (
            "IRED targets were already revealed by Phase 2D. "
            "These V9 results are development diagnostics only and cannot "
            "repair or replace the frozen Phase-2 verdict."
        ),
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "V9_PAIR_TRANSFER_RESULTS.json").write_text(
        json.dumps(result, indent=2, default=str),
        encoding="utf-8",
    )
    pd.DataFrame(
        {
            "candidate_id": test_ids,
            "mutation_count": [len(x) for x in test_mutations],
            "target": test_target,
            "strict_v83_scoreable": strict_scoreable,
            "all_main_supported": all_main,
            "PHASE2_V8_3_ABSTENTION": phase2_prediction,
            "B2_MAIN_ONLY": b2_prediction,
            "V9_PAIR_TRANSFER": v9_prediction,
            "V9_PAIR_TRANSFER_PERMUTED_CONTROL": permuted_prediction,
            "v9_exact_pair_count": diagnostics["exact_pair_count"],
            "v9_transferred_pair_count": diagnostics["transferred_pair_count"],
            "v9_unsupported_pair_count": diagnostics["unsupported_pair_count"],
        }
    ).to_csv(
        output_dir / "V9_PAIR_TRANSFER_PREDICTIONS.csv",
        index=False,
    )
    print(json.dumps(result, indent=2, default=str))
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source_gz")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    run(Path(args.source_gz), Path(args.out))


if __name__ == "__main__":
    main()
