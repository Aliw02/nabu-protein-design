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

from nabu_protein.higher_order import additive_score
from nabu_protein.v83 import NabuV83Model


SEED = 161
EXPECTED_IRED_SHA256 = (
    "aa45a2f85fb1af87b6b0e86397b3f8f292a061dc574be2536673b5bd490b0e74"
)
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
VERSION = "NABU_V9_0_FACTORIZED_PAIR_TRANSFER_DEV"


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
        raise ValueError("Fit sequences must have one fixed length.")
    if any(not set(sequence).issubset(STANDARD_AA) for sequence in sequences):
        raise ValueError("Fit sequences contain non-standard amino acids.")
    length = next(iter(lengths))
    return "".join(
        Counter(sequence[index] for sequence in sequences).most_common(1)[0][0]
        for index in range(length)
    )


def mutation_set(sequence: str, reference: str) -> tuple[str, ...]:
    sequence = str(sequence).strip().upper()
    if len(sequence) != len(reference):
        raise ValueError("Sequence/reference length mismatch.")
    if not set(sequence).issubset(STANDARD_AA):
        raise ValueError("Sequence contains non-standard amino acids.")
    return tuple(
        f"{source}{index + 1}{target}"
        for index, (source, target) in enumerate(zip(reference, sequence))
        if source != target
    )


def evaluate(target: np.ndarray, prediction: np.ndarray) -> dict:
    target = np.asarray(target, dtype=float)
    prediction = np.asarray(prediction, dtype=float)
    if len(target) != len(prediction) or len(target) == 0:
        raise ValueError("Target/prediction length mismatch.")
    if not np.isfinite(target).all() or not np.isfinite(prediction).all():
        raise ValueError("Metrics require finite values.")

    rho = float(spearmanr(target, prediction).statistic)
    relevance = target - float(np.min(target))
    ndcg = float(ndcg_score(relevance[None, :], prediction[None, :]))

    n = len(target)
    k = int(math.ceil(0.01 * n))
    predicted_top = np.argsort(-prediction, kind="stable")[:k]
    true_top = np.argsort(-target, kind="stable")[:k]
    hits = len(set(predicted_top.tolist()) & set(true_top.tolist()))

    global_min = float(np.min(target))
    global_max = float(np.max(target))
    best_selected = float(np.max(target[predicted_top]))
    span = global_max - global_min
    regret = 0.0 if span <= 0.0 else float(
        (global_max - best_selected) / span
    )

    return {
        "spearman": rho,
        "ndcg": ndcg,
        "top1_percent_k": k,
        "top1_percent_hits": int(hits),
        "top1_percent_recall": float(hits / k),
        "top1_percent_enrichment": float(hits * n / (k * k)),
        "normalized_regret_top1pct": regret,
        "best_true_in_predicted_top1pct": best_selected,
    }


def fit_factorized_pair_transfer(
    mutation_sets: list[tuple[str, ...]],
    labels: np.ndarray,
    oof_b2: np.ndarray,
    permute: bool = False,
) -> dict:
    labels = np.asarray(labels, dtype=float)
    oof_b2 = np.asarray(oof_b2, dtype=float)
    double_indices = np.array(
        [index for index, value in enumerate(mutation_sets) if len(value) == 2],
        dtype=int,
    )
    if len(double_indices) == 0:
        raise RuntimeError("No double-mutant rows available for pair transfer.")

    pair_rows = [mutation_sets[index] for index in double_indices]
    residual = labels[double_indices] - oof_b2[double_indices]

    nodes = sorted({mutation for pair in pair_rows for mutation in pair})
    node_to_col = {mutation: index + 1 for index, mutation in enumerate(nodes)}
    design = np.zeros((len(pair_rows), len(nodes) + 1), dtype=float)
    design[:, 0] = 1.0
    support = Counter()

    for row_index, pair in enumerate(pair_rows):
        for mutation in pair:
            design[row_index, node_to_col[mutation]] = 1.0
            support[mutation] += 1

    target = residual.copy()
    if permute:
        target = np.random.default_rng(SEED).permutation(target)

    coefficient, _, rank, singular_values = np.linalg.lstsq(
        design,
        target,
        rcond=None,
    )

    node_effect = {
        mutation: float(coefficient[node_to_col[mutation]])
        for mutation in nodes
    }
    node_confidence = {
        mutation: float(support[mutation] / (support[mutation] + 1.0))
        for mutation in nodes
    }

    fitted = design @ coefficient
    fit_rmse = float(np.sqrt(np.mean((target - fitted) ** 2)))

    return {
        "bias": float(coefficient[0]),
        "node_effect": node_effect,
        "node_confidence": node_confidence,
        "node_support": {key: int(value) for key, value in support.items()},
        "double_row_count": int(len(pair_rows)),
        "node_count": int(len(nodes)),
        "matrix_rank": int(rank),
        "singular_value_min": float(np.min(singular_values)),
        "singular_value_max": float(np.max(singular_values)),
        "fit_rmse": fit_rmse,
        "permuted": bool(permute),
    }


def transferred_pair_value(pair: tuple[str, str], factor_model: dict) -> float | None:
    left, right = pair
    effects = factor_model["node_effect"]
    confidence = factor_model["node_confidence"]
    if left not in effects or right not in effects:
        return None
    raw = (
        float(factor_model["bias"])
        + float(effects[left])
        + float(effects[right])
    )
    weight = math.sqrt(float(confidence[left]) * float(confidence[right]))
    return float(raw * weight)


def score_v9(
    mutation_sets: list[tuple[str, ...]],
    frozen_model: NabuV83Model,
    factor_model: dict,
) -> tuple[np.ndarray, list[dict]]:
    base = frozen_model.model["base"]
    global_mean = float(base["global_mean"])
    main_memory = base["main"]
    pair_memory = base["pair"]

    prediction = np.empty(len(mutation_sets), dtype=float)
    diagnostics = []

    for index, mutations in enumerate(mutation_sets):
        all_main = all(mutation in main_memory for mutation in mutations)
        if not all_main:
            prediction[index] = global_mean
            diagnostics.append(
                {
                    "all_main_supported": False,
                    "exact_pair_count": 0,
                    "transferred_pair_count": 0,
                    "missing_factor_pair_count": 0,
                }
            )
            continue

        score = float(additive_score(mutations, global_mean, main_memory))
        exact_count = 0
        transfer_count = 0
        missing_factor_count = 0

        for pair in combinations(mutations, 2):
            exact = pair_memory.get(pair)
            if exact is not None:
                score += float(exact["mean"]) * float(exact["confidence"])
                exact_count += 1
                continue

            transferred = transferred_pair_value(pair, factor_model)
            if transferred is None:
                missing_factor_count += 1
                continue
            score += transferred
            transfer_count += 1

        prediction[index] = score
        diagnostics.append(
            {
                "all_main_supported": True,
                "exact_pair_count": exact_count,
                "transferred_pair_count": transfer_count,
                "missing_factor_pair_count": missing_factor_count,
            }
        )

    return prediction, diagnostics


def score_b2(
    mutation_sets: list[tuple[str, ...]],
    frozen_model: NabuV83Model,
) -> tuple[np.ndarray, np.ndarray]:
    base = frozen_model.model["base"]
    global_mean = float(base["global_mean"])
    main_memory = base["main"]
    prediction = np.empty(len(mutation_sets), dtype=float)
    all_main = np.zeros(len(mutation_sets), dtype=bool)

    for index, mutations in enumerate(mutation_sets):
        supported = all(mutation in main_memory for mutation in mutations)
        all_main[index] = supported
        prediction[index] = (
            float(additive_score(mutations, global_mean, main_memory))
            if supported
            else global_mean
        )
    return prediction, all_main


def score_phase2_abstention(
    mutation_sets: list[tuple[str, ...]],
    candidate_ids: list[str],
    frozen_model: NabuV83Model,
) -> tuple[np.ndarray, np.ndarray]:
    base = frozen_model.model["base"]
    main_memory = base["main"]
    pair_memory = base["pair"]

    strict = np.array(
        [
            all(mutation in main_memory for mutation in mutations)
            and any(pair in pair_memory for pair in combinations(mutations, 2))
            for mutations in mutation_sets
        ],
        dtype=bool,
    )
    if int(strict.sum()) == 0:
        raise RuntimeError("No strict V8.3-scoreable rows.")

    strict_mutations = [
        mutations for mutations, use in zip(mutation_sets, strict) if use
    ]
    strict_ids = [
        candidate_id for candidate_id, use in zip(candidate_ids, strict) if use
    ]
    strict_frame = frozen_model.score_candidates(strict_mutations, strict_ids)
    strict_prediction = strict_frame[
        "V8_3_ADAPTIVE_ROUTER"
    ].to_numpy(dtype=float)
    if not np.isfinite(strict_prediction).all():
        raise RuntimeError("Frozen V8.3 produced non-finite strict scores.")

    floor = float(np.min(strict_prediction) - 1.0)
    prediction = np.full(len(mutation_sets), floor, dtype=float)
    prediction[strict] = strict_prediction
    return prediction, strict


def subset_spearman(
    target: np.ndarray,
    prediction: np.ndarray,
    mask: np.ndarray,
) -> float:
    if int(mask.sum()) < 2:
        return float("nan")
    return float(spearmanr(target[mask], prediction[mask]).statistic)


def per_mutation_count(
    mutation_sets: list[tuple[str, ...]],
    target: np.ndarray,
    arms: dict[str, np.ndarray],
) -> dict:
    counts = np.array([len(value) for value in mutation_sets], dtype=int)
    result = {}
    for mutation_count in sorted(set(counts.tolist())):
        mask = counts == mutation_count
        if int(mask.sum()) < 2:
            continue
        result[str(mutation_count)] = {
            "count": int(mask.sum()),
            "spearman": {
                name: float(spearmanr(target[mask], score[mask]).statistic)
                for name, score in arms.items()
            },
        }
    return result


def run(source_gz: Path, output_dir: Path) -> dict:
    source_sha = sha256_file(source_gz)
    if source_sha != EXPECTED_IRED_SHA256:
        raise RuntimeError(
            f"IRED SHA256 mismatch: expected {EXPECTED_IRED_SHA256}, got {source_sha}"
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

    fit = frame[set_text.eq("train") & ~validation_true].copy()
    validation = frame[validation_true].copy()
    test = frame[set_text.eq("test")].copy()
    if (len(fit), len(validation), len(test)) != (3746, 662, 4178):
        raise RuntimeError(
            f"Unexpected IRED split counts: {len(fit)}/{len(validation)}/{len(test)}"
        )

    reference = derive_reference(fit["sequence"].tolist())
    fit_mutations = [mutation_set(value, reference) for value in fit["sequence"]]
    test_mutations = [mutation_set(value, reference) for value in test["sequence"]]

    fit_labels = fit["target"].to_numpy(dtype=float)
    test_labels = test["target"].to_numpy(dtype=float)
    model = NabuV83Model().fit(
        fit_mutations,
        fit_labels,
        fit["sequence"].tolist(),
    )

    phase2_prediction, strict_mask = score_phase2_abstention(
        test_mutations,
        test["sequence"].tolist(),
        model,
    )
    b2_prediction, all_main_mask = score_b2(test_mutations, model)

    pair_transfer = fit_factorized_pair_transfer(
        fit_mutations,
        fit_labels,
        np.asarray(model.model["oof_b2"], dtype=float),
        permute=False,
    )
    pair_transfer_permuted = fit_factorized_pair_transfer(
        fit_mutations,
        fit_labels,
        np.asarray(model.model["oof_b2"], dtype=float),
        permute=True,
    )

    v9_prediction, v9_diag = score_v9(
        test_mutations,
        model,
        pair_transfer,
    )
    perm_prediction, _ = score_v9(
        test_mutations,
        model,
        pair_transfer_permuted,
    )

    # Deterministic replay of the new arm.
    replay_prediction, _ = score_v9(
        test_mutations,
        model,
        pair_transfer,
    )
    deterministic = bool(np.array_equal(v9_prediction, replay_prediction))

    arms = {
        "PHASE2_V8_3_ABSTENTION": phase2_prediction,
        "B2_MAIN_ONLY": b2_prediction,
        "V9_PAIR_TRANSFER": v9_prediction,
        "V9_PAIR_TRANSFER_PERMUTED_CONTROL": perm_prediction,
    }
    if not all(np.isfinite(value).all() for value in arms.values()):
        raise RuntimeError("At least one arm produced non-finite predictions.")

    metrics = {name: evaluate(test_labels, score) for name, score in arms.items()}

    known_main_no_exact_pair = all_main_mask & ~strict_mask
    non_strict_transfer = {
        "count": int(known_main_no_exact_pair.sum()),
        "B2_MAIN_ONLY_spearman": subset_spearman(
            test_labels, b2_prediction, known_main_no_exact_pair
        ),
        "V9_PAIR_TRANSFER_spearman": subset_spearman(
            test_labels, v9_prediction, known_main_no_exact_pair
        ),
    }
    non_strict_transfer["delta_spearman"] = float(
        non_strict_transfer["V9_PAIR_TRANSFER_spearman"]
        - non_strict_transfer["B2_MAIN_ONLY_spearman"]
    )

    v9_metric = metrics["V9_PAIR_TRANSFER"]
    b2_metric = metrics["B2_MAIN_ONLY"]
    perm_metric = metrics["V9_PAIR_TRANSFER_PERMUTED_CONTROL"]

    decision_checks = {
        "spearman_gt_b2": bool(
            v9_metric["spearman"] > b2_metric["spearman"]
        ),
        "elite_hits_gt_b2": bool(
            v9_metric["top1_percent_hits"] > b2_metric["top1_percent_hits"]
        ),
        "regret_not_worse_than_b2": bool(
            v9_metric["normalized_regret_top1pct"]
            <= b2_metric["normalized_regret_top1pct"]
        ),
        "spearman_gt_permuted_control": bool(
            v9_metric["spearman"] > perm_metric["spearman"]
        ),
        "finite_predictions": True,
        "deterministic_replay": deterministic,
        "gain_outside_original_strict_subset": bool(
            non_strict_transfer["delta_spearman"] > 0.0
        ),
    }
    promising = bool(
        decision_checks["spearman_gt_b2"]
        and decision_checks["elite_hits_gt_b2"]
        and decision_checks["regret_not_worse_than_b2"]
        and decision_checks["spearman_gt_permuted_control"]
        and decision_checks["finite_predictions"]
        and decision_checks["deterministic_replay"]
        and decision_checks["gain_outside_original_strict_subset"]
    )

    diag_frame = pd.DataFrame(v9_diag)
    counts = np.array([len(value) for value in test_mutations], dtype=int)

    output_dir.mkdir(parents=True, exist_ok=True)
    predictions = pd.DataFrame(
        {
            "candidate_id": test["sequence"].tolist(),
            "mutation_count": counts,
            "target": test_labels,
            "strict_v83_scoreable": strict_mask,
            "all_main_supported": all_main_mask,
            **arms,
            "v9_exact_pair_count": diag_frame["exact_pair_count"].to_numpy(),
            "v9_transferred_pair_count": diag_frame[
                "transferred_pair_count"
            ].to_numpy(),
            "v9_missing_factor_pair_count": diag_frame[
                "missing_factor_pair_count"
            ].to_numpy(),
        }
    )
    predictions.to_csv(output_dir / "V9_IRED_PREDICTIONS.csv", index=False)

    result = {
        "version": VERSION,
        "status": (
            "MECHANISTICALLY_PROMISING" if promising else "DEVELOPMENT_FAIL"
        ),
        "scientific_claim": False,
        "phase3_opened": False,
        "source_sha256": source_sha,
        "frozen_v83_unchanged": True,
        "fit_count": int(len(fit)),
        "validation_held_out_count": int(len(validation)),
        "test_diagnostic_count": int(len(test)),
        "reference_length": int(len(reference)),
        "coverage": {
            "all_main_supported": int(all_main_mask.sum()),
            "strict_v83_scoreable": int(strict_mask.sum()),
            "known_main_without_exact_pair_scoreability": int(
                known_main_no_exact_pair.sum()
            ),
        },
        "factor_model": {
            key: value
            for key, value in pair_transfer.items()
            if key not in {"node_effect", "node_confidence", "node_support"}
        },
        "permuted_factor_model": {
            key: value
            for key, value in pair_transfer_permuted.items()
            if key not in {"node_effect", "node_confidence", "node_support"}
        },
        "v9_pair_diagnostics": {
            "rows_with_any_transferred_pair": int(
                (diag_frame["transferred_pair_count"] > 0).sum()
            ),
            "total_exact_pair_contributions": int(
                diag_frame["exact_pair_count"].sum()
            ),
            "total_transferred_pair_contributions": int(
                diag_frame["transferred_pair_count"].sum()
            ),
            "total_missing_factor_pairs": int(
                diag_frame["missing_factor_pair_count"].sum()
            ),
            "factor_node_count": int(pair_transfer["node_count"]),
        },
        "metrics": metrics,
        "known_main_no_exact_pair_diagnostic": non_strict_transfer,
        "per_mutation_count": per_mutation_count(
            test_mutations,
            test_labels,
            arms,
        ),
        "decision_checks": decision_checks,
        "decision": (
            "CONTINUE_V9_PAIR_TRANSFER"
            if promising
            else "DO_NOT_PROMOTE_V9_0"
        ),
        "interpretation_boundary": (
            "IRED labels were already revealed by Phase 2D. This V9 result is "
            "development evidence only and cannot repair or replace the frozen "
            "Phase-2 scientific verdict. Phase 3 remains unopened."
        ),
    }

    (output_dir / "V9_RESULTS.json").write_text(
        json.dumps(result, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2))
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("ired_gz")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    run(Path(args.ired_gz), Path(args.out))


if __name__ == "__main__":
    main()
