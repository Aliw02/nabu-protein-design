from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import ndcg_score

from nabu_protein.higher_order import FROZEN_CORE_COMMIT, model_diagnostics
from nabu_protein.v83 import NabuV83Model
from ired_scoreability_preflight import derive_reference, mutation_set, scoreability_preflight


EXPECTED_IRED_SHA256 = "aa45a2f85fb1af87b6b0e86397b3f8f292a061dc574be2536673b5bd490b0e74"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_predictions(frame: pd.DataFrame) -> str:
    payload = "\n".join(
        f"{row.candidate_id},{float(row.prediction):.17g}"
        for row in frame[["candidate_id", "prediction"]].itertuples(index=False)
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def evaluate_complete_test(target: np.ndarray, prediction: np.ndarray) -> dict:
    if len(target) != len(prediction) or len(target) == 0:
        raise ValueError("Target/prediction length mismatch or empty test.")
    if not np.isfinite(target).all() or not np.isfinite(prediction).all():
        raise ValueError("Metrics require finite full-test values.")

    rho = float(spearmanr(target, prediction).statistic)
    relevance = target - float(np.min(target))
    ndcg = float(ndcg_score(relevance[None, :], prediction[None, :]))

    n = len(target)
    k = int(math.ceil(0.01 * n))
    predicted_top = np.argsort(-prediction, kind="stable")[:k]
    true_top = np.argsort(-target, kind="stable")[:k]
    hits = len(set(predicted_top.tolist()).intersection(true_top.tolist()))

    top1_recall = float(hits / k)
    top1_enrichment = float(hits * n / (k * k))

    best_true_in_predicted_top = float(np.max(target[predicted_top]))
    global_min = float(np.min(target))
    global_max = float(np.max(target))
    span = global_max - global_min
    normalized_regret = (
        0.0
        if span <= 0.0
        else float((global_max - best_true_in_predicted_top) / span)
    )

    return {
        "spearman": rho,
        "ndcg": ndcg,
        "top1_percent_k": int(k),
        "top1_percent_hits": int(hits),
        "top1_percent_recall": top1_recall,
        "top1_percent_enrichment": top1_enrichment,
        "normalized_regret_top1pct": normalized_regret,
        "best_true_in_predicted_top1pct": best_true_in_predicted_top,
        "global_test_min": global_min,
        "global_test_max": global_max,
    }


def run(source_gz: Path, output_dir: Path) -> dict:
    actual_sha = sha256_file(source_gz)
    if actual_sha != EXPECTED_IRED_SHA256:
        raise RuntimeError(
            f"IRED source SHA256 mismatch: expected {EXPECTED_IRED_SHA256}, got {actual_sha}"
        )

    # Mandatory identity-only gate before this function reads any target value.
    identity_gate = scoreability_preflight(source_gz)
    if not identity_gate["full_test_scoreable"]:
        raise RuntimeError(
            "Frozen V8.3 cannot score the complete official IRED test set; "
            "target reveal is prohibited until a versioned full-test policy exists."
        )

    frame = pd.read_csv(source_gz, compression="gzip")
    required = {"sequence", "target", "set", "validation"}
    missing = required - set(frame.columns)
    if missing:
        raise RuntimeError(f"IRED source missing columns: {sorted(missing)}")

    frame["sequence"] = frame["sequence"].astype(str).str.strip().str.upper()
    frame["target"] = pd.to_numeric(frame["target"], errors="raise")
    if not np.isfinite(frame["target"].to_numpy(dtype=float)).all():
        raise RuntimeError("IRED contains non-finite target values.")

    set_text = frame["set"].astype(str).str.lower()
    validation_true = frame["validation"].fillna(False).astype(bool)
    fit_frame = frame[set_text.eq("train") & ~validation_true].copy()
    validation_frame = frame[validation_true].copy()
    test_frame = frame[set_text.eq("test")].copy()

    if (len(fit_frame), len(validation_frame), len(test_frame)) != (3746, 662, 4178):
        raise RuntimeError(
            "Unexpected first-touch IRED split counts: "
            f"{len(fit_frame)}/{len(validation_frame)}/{len(test_frame)}"
        )

    reference = derive_reference(fit_frame["sequence"].tolist())
    fit_mutations = [
        mutation_set(sequence, reference)
        for sequence in fit_frame["sequence"]
    ]
    test_mutations = [
        mutation_set(sequence, reference)
        for sequence in test_frame["sequence"]
    ]

    model = NabuV83Model().fit(
        mutation_sets=fit_mutations,
        labels=fit_frame["target"].to_numpy(dtype=float),
        candidate_ids=fit_frame["sequence"].tolist(),
    )
    scored = model.score_candidates(
        mutation_sets=test_mutations,
        candidate_ids=test_frame["sequence"].tolist(),
    )

    if not bool(scored["scoreable"].all()):
        raise RuntimeError(
            "Runtime scoreability disagrees with identity-only preflight."
        )
    prediction = scored["V8_3_ADAPTIVE_ROUTER"].to_numpy(dtype=float)
    target = test_frame["target"].to_numpy(dtype=float)
    if not np.isfinite(prediction).all():
        raise RuntimeError("Frozen NABU produced non-finite full-test predictions.")

    metrics = evaluate_complete_test(target, prediction)

    output_dir.mkdir(parents=True, exist_ok=True)
    predictions = pd.DataFrame(
        {
            "candidate_id": test_frame["sequence"].tolist(),
            "mutation_count": [len(value) for value in test_mutations],
            "target": target,
            "prediction": prediction,
        }
    )
    predictions["true_rank"] = (
        predictions["target"]
        .rank(method="min", ascending=False)
        .astype(int)
    )
    predictions["predicted_rank"] = (
        predictions["prediction"]
        .rank(method="min", ascending=False)
        .astype(int)
    )
    predictions.to_csv(output_dir / "IRED_TEST_PREDICTIONS.csv", index=False)

    manifest = {
        "version": "NABU_PHASE2D_IRED_BLIND_V1",
        "scientific_claim": False,
        "status": "BLIND_RESULT_COMPLETE",
        "source_sha256": actual_sha,
        "source_expected_sha256": EXPECTED_IRED_SHA256,
        "frozen_phase2c_version": "NABU_PHASE2C_MULTILANDSCAPE_V3",
        "frozen_core_commit": FROZEN_CORE_COMMIT,
        "fit_rule": "set == train AND validation != True",
        "validation_rule": "validation == True held out; no tuning",
        "test_rule": "set == test",
        "fit_count": int(len(fit_frame)),
        "validation_count": int(len(validation_frame)),
        "test_count": int(len(test_frame)),
        "reference_sha256": hashlib.sha256(reference.encode("utf-8")).hexdigest(),
        "reference_length": int(len(reference)),
        "identity_preflight": identity_gate,
        "router_decision": model.router_decision,
        "model_diagnostics": model_diagnostics(model.model),
        "prediction_sha256": sha256_predictions(predictions),
        "metrics": metrics,
        "primary_metrics": {
            "spearman": metrics["spearman"],
            "ndcg": metrics["ndcg"],
        },
        "secondary_metrics": {
            key: metrics[key]
            for key in (
                "top1_percent_k",
                "top1_percent_hits",
                "top1_percent_recall",
                "top1_percent_enrichment",
                "normalized_regret_top1pct",
            )
        },
        "interpretation_rule": (
            "Report this untouched result as observed. No Phase-2 scientific "
            "retuning is allowed after reveal."
        ),
    }
    (output_dir / "IRED_BLIND_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, default=str),
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2, default=str))
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source_gz")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    run(Path(args.source_gz), Path(args.out))


if __name__ == "__main__":
    main()
