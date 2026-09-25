from __future__ import annotations

import argparse
import hashlib
import json
import math
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import ndcg_score

from nabu_protein.higher_order import (
    FROZEN_CORE_COMMIT,
    additive_score,
    model_diagnostics,
)
from nabu_protein.v83 import NabuV83Model
from ired_scoreability_preflight import (
    derive_reference,
    mutation_set,
    scoreability_preflight,
)


EXPECTED_IRED_SHA256 = (
    "aa45a2f85fb1af87b6b0e86397b3f8f292a061dc574be2536673b5bd490b0e74"
)
BACKOFF_VERSION = "IRED_EVIDENCE_BACKOFF_V1"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_predictions(frame: pd.DataFrame) -> str:
    payload = "\n".join(
        f"{row.candidate_id},{float(row.prediction):.17g},{row.prediction_source}"
        for row in frame[
            ["candidate_id", "prediction", "prediction_source"]
        ].itertuples(index=False)
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def evaluate_complete_test(
    target: np.ndarray,
    prediction: np.ndarray,
) -> dict:
    if len(target) != len(prediction) or len(target) == 0:
        raise ValueError("Target/prediction length mismatch or empty test.")
    if not np.isfinite(target).all() or not np.isfinite(prediction).all():
        raise ValueError("Metrics require finite full-test values.")

    rho = float(spearmanr(target, prediction).statistic)
    relevance = target - float(np.min(target))
    ndcg = float(
        ndcg_score(
            relevance[None, :],
            prediction[None, :],
        )
    )

    n = len(target)
    k = int(math.ceil(0.01 * n))
    predicted_top = np.argsort(
        -prediction,
        kind="stable",
    )[:k]
    true_top = np.argsort(
        -target,
        kind="stable",
    )[:k]
    hits = len(
        set(predicted_top.tolist()).intersection(
            true_top.tolist()
        )
    )

    top1_recall = float(hits / k)
    top1_enrichment = float(hits * n / (k * k))

    best_true_in_predicted_top = float(
        np.max(target[predicted_top])
    )
    global_min = float(np.min(target))
    global_max = float(np.max(target))
    span = global_max - global_min
    normalized_regret = (
        0.0
        if span <= 0.0
        else float(
            (global_max - best_true_in_predicted_top)
            / span
        )
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


def predict_with_frozen_backoff(
    model: NabuV83Model,
    mutation_sets: list[tuple[str, ...]],
    candidate_ids: list[str],
) -> pd.DataFrame:
    if model.model is None:
        raise RuntimeError("Model must be fit before prediction.")
    if len(mutation_sets) != len(candidate_ids):
        raise ValueError(
            "mutation_sets and candidate_ids must have equal length."
        )

    base = model.model["base"]
    global_mean = float(base["global_mean"])
    main_memory = base["main"]
    pair_memory = base["pair"]

    rows = []
    strict_mutations = []
    strict_ids = []
    strict_positions = []

    for index, (mutations, candidate_id) in enumerate(
        zip(mutation_sets, candidate_ids)
    ):
        all_main_supported = all(
            mutation in main_memory
            for mutation in mutations
        )
        supported_pair_count = sum(
            pair in pair_memory
            for pair in combinations(mutations, 2)
        )
        any_pair_supported = supported_pair_count > 0

        if all_main_supported and any_pair_supported:
            prediction_source = "V8_3_ADAPTIVE_ROUTER"
            prediction = np.nan
            strict_positions.append(index)
            strict_mutations.append(mutations)
            strict_ids.append(candidate_id)
        elif all_main_supported:
            prediction_source = "B2_ADDITIVE"
            prediction = float(
                additive_score(
                    mutations,
                    global_mean,
                    main_memory,
                )
            )
        else:
            prediction_source = "GLOBAL_MEAN_ABSTENTION"
            prediction = global_mean

        rows.append(
            {
                "candidate_id": str(candidate_id),
                "all_main_supported": bool(
                    all_main_supported
                ),
                "supported_pair_count": int(
                    supported_pair_count
                ),
                "prediction_source": prediction_source,
                "prediction": prediction,
            }
        )

    if strict_ids:
        strict_scored = model.score_candidates(
            strict_mutations,
            strict_ids,
        )
        if not bool(strict_scored["scoreable"].all()):
            raise RuntimeError(
                "Frozen router path disagrees with structural "
                "scoreability classification."
            )
        strict_prediction = strict_scored[
            "V8_3_ADAPTIVE_ROUTER"
        ].to_numpy(dtype=float)
        if not np.isfinite(strict_prediction).all():
            raise RuntimeError(
                "Frozen router produced non-finite strict predictions."
            )
        for index, value in zip(
            strict_positions,
            strict_prediction,
        ):
            rows[index]["prediction"] = float(value)

    result = pd.DataFrame(rows)
    if not np.isfinite(
        result["prediction"].to_numpy(dtype=float)
    ).all():
        raise RuntimeError(
            "Evidence backoff produced non-finite prediction."
        )
    return result


def run(source_gz: Path, output_dir: Path) -> dict:
    actual_sha = sha256_file(source_gz)
    if actual_sha != EXPECTED_IRED_SHA256:
        raise RuntimeError(
            "IRED source SHA256 mismatch: "
            f"expected {EXPECTED_IRED_SHA256}, got {actual_sha}"
        )

    # Mandatory label-blind identity gate. This reads no target values.
    identity_gate = scoreability_preflight(source_gz)
    if (
        identity_gate["fit_count"],
        identity_gate["validation_count"],
        identity_gate["test_count"],
    ) != (3746, 662, 4178):
        raise RuntimeError(
            "Identity preflight split counts differ from frozen contract."
        )

    # Target reveal starts only after all structural checks above.
    frame = pd.read_csv(source_gz, compression="gzip")
    required = {
        "sequence",
        "target",
        "set",
        "validation",
    }
    missing = required - set(frame.columns)
    if missing:
        raise RuntimeError(
            f"IRED source missing columns: {sorted(missing)}"
        )

    frame["sequence"] = (
        frame["sequence"]
        .astype(str)
        .str.strip()
        .str.upper()
    )
    frame["target"] = pd.to_numeric(
        frame["target"],
        errors="raise",
    )
    if not np.isfinite(
        frame["target"].to_numpy(dtype=float)
    ).all():
        raise RuntimeError(
            "IRED contains non-finite target values."
        )

    set_text = frame["set"].astype(str).str.lower()
    validation_true = (
        frame["validation"]
        .fillna(False)
        .astype(bool)
    )
    fit_frame = frame[
        set_text.eq("train") & ~validation_true
    ].copy()
    validation_frame = frame[
        validation_true
    ].copy()
    test_frame = frame[
        set_text.eq("test")
    ].copy()

    if (
        len(fit_frame),
        len(validation_frame),
        len(test_frame),
    ) != (3746, 662, 4178):
        raise RuntimeError(
            "Unexpected frozen IRED split counts: "
            f"{len(fit_frame)}/"
            f"{len(validation_frame)}/"
            f"{len(test_frame)}"
        )

    reference = derive_reference(
        fit_frame["sequence"].tolist()
    )
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
        labels=fit_frame[
            "target"
        ].to_numpy(dtype=float),
        candidate_ids=fit_frame[
            "sequence"
        ].tolist(),
    )

    predicted = predict_with_frozen_backoff(
        model=model,
        mutation_sets=test_mutations,
        candidate_ids=test_frame[
            "sequence"
        ].tolist(),
    )

    source_counts = (
        predicted["prediction_source"]
        .value_counts()
        .to_dict()
    )
    if (
        int(
            source_counts.get(
                "V8_3_ADAPTIVE_ROUTER",
                0,
            )
        )
        != identity_gate["scoreable_count"]
    ):
        raise RuntimeError(
            "Runtime strict-router count disagrees "
            "with label-blind preflight."
        )

    prediction = predicted[
        "prediction"
    ].to_numpy(dtype=float)
    target = test_frame[
        "target"
    ].to_numpy(dtype=float)

    metrics = evaluate_complete_test(
        target,
        prediction,
    )

    strict_mask = (
        predicted["prediction_source"]
        == "V8_3_ADAPTIVE_ROUTER"
    ).to_numpy()
    strict_metrics = None
    if int(strict_mask.sum()) >= 2:
        strict_target = target[strict_mask]
        strict_prediction = prediction[
            strict_mask
        ]
        strict_metrics = {
            "count": int(strict_mask.sum()),
            "spearman": float(
                spearmanr(
                    strict_target,
                    strict_prediction,
                ).statistic
            ),
            "ndcg": float(
                ndcg_score(
                    (
                        strict_target
                        - float(
                            np.min(
                                strict_target
                            )
                        )
                    )[None, :],
                    strict_prediction[
                        None,
                        :
                    ],
                )
            ),
            "role": "DIAGNOSTIC_ONLY",
        }

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )
    predictions = pd.DataFrame(
        {
            "candidate_id": test_frame[
                "sequence"
            ].tolist(),
            "mutation_count": [
                len(value)
                for value in test_mutations
            ],
            "target": target,
            "prediction": prediction,
            "prediction_source": predicted[
                "prediction_source"
            ].tolist(),
            "all_main_supported": predicted[
                "all_main_supported"
            ].tolist(),
            "supported_pair_count": predicted[
                "supported_pair_count"
            ].tolist(),
        }
    )
    predictions["true_rank"] = (
        predictions["target"]
        .rank(
            method="min",
            ascending=False,
        )
        .astype(int)
    )
    predictions["predicted_rank"] = (
        predictions["prediction"]
        .rank(
            method="min",
            ascending=False,
        )
        .astype(int)
    )
    predictions.to_csv(
        output_dir
        / "IRED_TEST_PREDICTIONS.csv",
        index=False,
    )

    manifest = {
        "version": (
            "NABU_PHASE2D_IRED_BLIND_V2"
        ),
        "scientific_claim": False,
        "status": "BLIND_RESULT_COMPLETE",
        "source_sha256": actual_sha,
        "source_expected_sha256": (
            EXPECTED_IRED_SHA256
        ),
        "frozen_phase2c_version": (
            "NABU_PHASE2C_"
            "MULTILANDSCAPE_V3"
        ),
        "frozen_core_commit": (
            FROZEN_CORE_COMMIT
        ),
        "full_test_policy": BACKOFF_VERSION,
        "fit_rule": (
            "set == train AND "
            "validation != True"
        ),
        "validation_rule": (
            "validation == True "
            "held out; no tuning"
        ),
        "test_rule": "set == test",
        "fit_count": int(
            len(fit_frame)
        ),
        "validation_count": int(
            len(validation_frame)
        ),
        "test_count": int(
            len(test_frame)
        ),
        "reference_sha256": (
            hashlib.sha256(
                reference.encode("utf-8")
            ).hexdigest()
        ),
        "reference_length": int(
            len(reference)
        ),
        "identity_preflight": identity_gate,
        "prediction_source_counts": {
            str(key): int(value)
            for key, value
            in source_counts.items()
        },
        "router_decision": (
            model.router_decision
        ),
        "model_diagnostics": (
            model_diagnostics(
                model.model
            )
        ),
        "prediction_sha256": (
            sha256_predictions(
                predictions
            )
        ),
        "metrics": metrics,
        "primary_metrics": {
            "spearman": (
                metrics["spearman"]
            ),
            "ndcg": (
                metrics["ndcg"]
            ),
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
        "strict_scoreable_diagnostic": (
            strict_metrics
        ),
        "interpretation_rule": (
            "Primary metrics use all "
            "4,178 official test rows "
            "under the preregistered "
            "identity-only evidence "
            "backoff. The strict "
            "scoreable subset is "
            "diagnostic only. No "
            "Phase-2 scientific "
            "retuning is allowed "
            "after reveal."
        ),
    }
    (
        output_dir
        / "IRED_BLIND_MANIFEST.json"
    ).write_text(
        json.dumps(
            manifest,
            indent=2,
            default=str,
        ),
        encoding="utf-8",
    )
    print(
        json.dumps(
            manifest,
            indent=2,
            default=str,
        )
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source_gz")
    parser.add_argument(
        "--out",
        required=True,
    )
    args = parser.parse_args()
    run(
        Path(args.source_gz),
        Path(args.out),
    )


if __name__ == "__main__":
    main()
