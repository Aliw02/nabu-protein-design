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

from nabu_protein.higher_order import score_base
from nabu_protein.v83 import NabuV83Model


EXPECTED_IRED_SHA256 = (
    "aa45a2f85fb1af87b6b0e86397b3f8f292a061dc574be2536673b5bd490b0e74"
)
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
VERSION = "NABU_V9_1_EVIDENCE_DEPTH_ROUTER_DEV_V1"


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
    low = float(np.min(target))
    high = float(np.max(target))
    best_true = float(np.max(target[predicted_top]))
    span = high - low

    return {
        "count": int(n),
        "spearman": finite_spearman(target, prediction),
        "ndcg": float(ndcg_score(relevance[None, :], prediction[None, :])),
        "top1_percent_k": int(k),
        "top1_percent_hits": hits,
        "top1_percent_recall": float(hits / k),
        "top1_percent_enrichment": float(hits * n / (k * k)),
        "normalized_regret_top1pct": (
            0.0 if span <= 0.0 else float((high - best_true) / span)
        ),
        "best_true_in_predicted_top1pct": best_true,
    }


def prediction_hash(candidate_ids, prediction) -> str:
    payload = "\n".join(
        f"{candidate_id},{float(value):.17g}"
        for candidate_id, value in zip(candidate_ids, prediction)
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def decide_pair_depth(model: NabuV83Model, fit_target) -> dict:
    target = np.asarray(fit_target, dtype=float)
    oof_b2 = np.asarray(model.model["oof_b2"], dtype=float)
    oof_b3 = np.asarray(model.model["oof_b3"], dtype=float)
    rho_b2 = float(spearmanr(oof_b2, target).statistic)
    rho_b3 = float(spearmanr(oof_b3, target).statistic)
    return {
        "oof_b2_spearman": rho_b2,
        "oof_b3_spearman": rho_b3,
        "b3_minus_b2": float(rho_b3 - rho_b2),
        "pair_depth_enabled": bool(rho_b3 > rho_b2),
        "rule": "enable exact B3 pair depth iff OOF_B3 Spearman > OOF_B2 Spearman",
    }


def b2_b3_predictions(model: NabuV83Model, mutation_sets):
    base = model.model["base"]
    main_memory = base["main"]
    pair_memory = base["pair"]
    global_mean = float(base["global_mean"])

    n = len(mutation_sets)
    b2 = np.full(n, global_mean, dtype=float)
    b3 = np.full(n, global_mean, dtype=float)
    all_main = np.zeros(n, dtype=bool)
    exact_pair = np.zeros(n, dtype=bool)
    exact_pair_count = np.zeros(n, dtype=int)

    for index, mutations in enumerate(mutation_sets):
        supported = all(mutation in main_memory for mutation in mutations)
        all_main[index] = supported
        if not supported:
            continue

        row_b2, row_b3 = score_base(mutations, base)
        b2[index] = float(row_b2)
        b3[index] = float(row_b3)

        pair_count = sum(
            pair in pair_memory
            for pair in combinations(mutations, 2)
        )
        exact_pair_count[index] = int(pair_count)
        exact_pair[index] = pair_count > 0

    return {
        "b2": b2,
        "b3": b3,
        "all_main": all_main,
        "exact_pair": exact_pair,
        "exact_pair_count": exact_pair_count,
    }


def v91_predictions(model, mutation_sets, fit_target):
    evidence = b2_b3_predictions(model, mutation_sets)
    gate = decide_pair_depth(model, fit_target)
    prediction = (
        evidence["b3"].copy()
        if gate["pair_depth_enabled"]
        else evidence["b2"].copy()
    )
    return prediction, evidence, gate


def phase2_abstention_predictions(model, mutation_sets, candidate_ids):
    base = model.model["base"]
    main_memory = base["main"]
    pair_memory = base["pair"]

    strict = np.array(
        [
            all(mutation in main_memory for mutation in mutations)
            and any(
                pair in pair_memory
                for pair in combinations(mutations, 2)
            )
            for mutations in mutation_sets
        ],
        dtype=bool,
    )
    if int(strict.sum()) == 0:
        raise RuntimeError("No strict V8.3-scoreable rows.")

    strict_indices = np.where(strict)[0]
    strict_scored = model.score_candidates(
        [mutation_sets[index] for index in strict_indices],
        [candidate_ids[index] for index in strict_indices],
    )
    values = strict_scored["V8_3_ADAPTIVE_ROUTER"].to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise RuntimeError("Frozen V8.3 strict score contains non-finite values.")

    floor = float(np.min(values) - 1.0)
    prediction = np.full(len(mutation_sets), floor, dtype=float)
    prediction[strict_indices] = values
    return prediction, strict


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


def masked_metrics(target, arms, mask):
    if int(mask.sum()) < 2:
        return None
    return {
        "count": int(mask.sum()),
        **{
            name: evaluate(target[mask], prediction[mask])
            for name, prediction in arms.items()
        },
    }


def run(source_gz: Path, output_dir: Path) -> dict:
    source_sha, fit_frame, validation_frame, test_frame = load_ired(source_gz)
    reference = derive_reference(fit_frame["sequence"].tolist())

    fit_mutations = [mutation_set(value, reference) for value in fit_frame["sequence"]]
    validation_mutations = [
        mutation_set(value, reference) for value in validation_frame["sequence"]
    ]
    test_mutations = [
        mutation_set(value, reference) for value in test_frame["sequence"]
    ]

    fit_target = fit_frame["target"].to_numpy(dtype=float)
    validation_target = validation_frame["target"].to_numpy(dtype=float)
    test_target = test_frame["target"].to_numpy(dtype=float)

    model = NabuV83Model().fit(
        mutation_sets=fit_mutations,
        labels=fit_target,
        candidate_ids=fit_frame["sequence"].tolist(),
    )

    test_ids = test_frame["sequence"].tolist()
    phase2_prediction, strict_mask = phase2_abstention_predictions(
        model,
        test_mutations,
        test_ids,
    )
    v91_prediction, test_evidence, gate = v91_predictions(
        model,
        test_mutations,
        fit_target,
    )

    b2_prediction = test_evidence["b2"]
    b3_prediction = test_evidence["b3"]
    all_main = test_evidence["all_main"]
    exact_pair = test_evidence["exact_pair"]

    replay_prediction, _, replay_gate = v91_predictions(
        model,
        test_mutations,
        fit_target,
    )
    deterministic = bool(
        np.array_equal(v91_prediction, replay_prediction)
        and gate == replay_gate
    )

    arms = {
        "PHASE2_V8_3_ABSTENTION": phase2_prediction,
        "B2_MAIN_ONLY": b2_prediction,
        "B3_AVAILABLE_ALWAYS": b3_prediction,
        "V9_1_EVIDENCE_DEPTH_ROUTER": v91_prediction,
    }
    if not all(np.isfinite(value).all() for value in arms.values()):
        raise RuntimeError("At least one V9.1 arm produced non-finite values.")

    metrics = {
        name: evaluate(test_target, prediction)
        for name, prediction in arms.items()
    }

    strict_all_main = all_main & exact_pair
    non_strict_all_main = all_main & ~exact_pair

    strict_metrics = masked_metrics(
        test_target,
        {
            "B2_MAIN_ONLY": b2_prediction,
            "B3_AVAILABLE_ALWAYS": b3_prediction,
            "V9_1_EVIDENCE_DEPTH_ROUTER": v91_prediction,
        },
        strict_all_main,
    )
    non_strict_metrics = masked_metrics(
        test_target,
        {
            "B2_MAIN_ONLY": b2_prediction,
            "B3_AVAILABLE_ALWAYS": b3_prediction,
            "V9_1_EVIDENCE_DEPTH_ROUTER": v91_prediction,
        },
        non_strict_all_main,
    )

    validation_v91, validation_evidence, validation_gate = v91_predictions(
        model,
        validation_mutations,
        fit_target,
    )
    validation_metrics = {
        "all_main_supported": int(validation_evidence["all_main"].sum()),
        "exact_pair_supported": int(validation_evidence["exact_pair"].sum()),
        "gate_matches_test": bool(validation_gate == gate),
        "B2_MAIN_ONLY": evaluate(
            validation_target,
            validation_evidence["b2"],
        ),
        "B3_AVAILABLE_ALWAYS": evaluate(
            validation_target,
            validation_evidence["b3"],
        ),
        "V9_1_EVIDENCE_DEPTH_ROUTER": evaluate(
            validation_target,
            validation_v91,
        ),
    }

    b2_metric = metrics["B2_MAIN_ONLY"]
    v91_metric = metrics["V9_1_EVIDENCE_DEPTH_ROUTER"]
    checks = {
        "finite_predictions": True,
        "deterministic_replay": deterministic,
        "coverage_expands_beyond_strict_v83": bool(
            int(all_main.sum()) > int(strict_mask.sum())
        ),
        "spearman_not_worse_than_b2": bool(
            v91_metric["spearman"] >= b2_metric["spearman"]
        ),
        "top1_hits_not_worse_than_b2": bool(
            v91_metric["top1_percent_hits"]
            >= b2_metric["top1_percent_hits"]
        ),
        "normalized_regret_not_worse_than_b2": bool(
            v91_metric["normalized_regret_top1pct"]
            <= b2_metric["normalized_regret_top1pct"]
        ),
    }
    retain = bool(all(checks.values()))

    output_dir.mkdir(parents=True, exist_ok=True)
    predictions = pd.DataFrame(
        {
            "candidate_id": test_ids,
            "mutation_count": [len(value) for value in test_mutations],
            "target": test_target,
            "strict_v83_scoreable": strict_mask,
            "all_main_supported": all_main,
            "exact_pair_supported": exact_pair,
            "exact_pair_count": test_evidence["exact_pair_count"],
            **arms,
        }
    )
    predictions.to_csv(
        output_dir / "V9_1_PREDICTIONS.csv",
        index=False,
    )

    result = {
        "version": VERSION,
        "status": "DEVELOPMENT_DIAGNOSTIC_AFTER_PHASE2_REVEAL",
        "scientific_claim": False,
        "phase3_opened": False,
        "source_sha256": source_sha,
        "frozen_v83_unchanged": True,
        "split_counts": {
            "fit": int(len(fit_frame)),
            "validation": int(len(validation_frame)),
            "test": int(len(test_frame)),
        },
        "reference_sha256": hashlib.sha256(
            reference.encode("utf-8")
        ).hexdigest(),
        "reference_length": int(len(reference)),
        "pair_depth_gate": gate,
        "coverage": {
            "strict_v83_scoreable": int(strict_mask.sum()),
            "all_main_supported": int(all_main.sum()),
            "exact_pair_supported_all_main": int(strict_all_main.sum()),
            "all_main_without_exact_pair": int(non_strict_all_main.sum()),
        },
        "metrics": metrics,
        "strict_exact_pair_diagnostic": strict_metrics,
        "non_strict_all_main_diagnostic": non_strict_metrics,
        "validation_diagnostic": validation_metrics,
        "prediction_sha256": {
            name: prediction_hash(test_ids, prediction)
            for name, prediction in arms.items()
        },
        "replay_sha256": prediction_hash(
            test_ids,
            replay_prediction,
        ),
        "decision_checks": checks,
        "development_decision": (
            "RETAIN_V9_1_AS_NEXT_ARCHITECTURE_BASE"
            if retain
            else "DO_NOT_PROMOTE_V9_1"
        ),
        "interpretation_boundary": (
            "IRED was revealed in Phase 2D and is development-only here. "
            "V9.1 changes scoreability/backoff semantics but does not alter "
            "frozen V8.3 or the Phase-2 verdict. Phase 3 remains unopened."
        ),
    }
    (output_dir / "V9_1_RESULTS.json").write_text(
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
