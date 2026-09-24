import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import pandas as pd


BUDGETS = [5, 10, 20, 40]
HERE = Path(__file__).resolve()
ROOT = HERE.parents[2]

PASSIVE_FINAL_PATH = (
    ROOT
    / "nabu_v8_3_rhla_sample_efficiency_final"
    / "FINAL_RESULTS.json"
)
PASSIVE_STAGE_B_HELPER = (
    ROOT
    / "experiments"
    / "v8_3_rhla_sample_efficiency"
    / "02_stage_b_reveal_evaluate.py"
)

spec = importlib.util.spec_from_file_location(
    "rhla_passive_eval",
    PASSIVE_STAGE_B_HELPER,
)
if spec is None or spec.loader is None:
    raise RuntimeError(
        f"Could not import passive evaluator from {PASSIVE_STAGE_B_HELPER}"
    )
passive_eval = importlib.util.module_from_spec(spec)
spec.loader.exec_module(passive_eval)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_frozen_rankings(scores, freeze):
    for arm, group in freeze.groupby("arm"):
        expected = (
            scores.sort_values(
                [arm, "candidate_id"],
                ascending=[False, True],
            )
            .head(50)["candidate_id"]
            .astype(str)
            .tolist()
        )
        actual = (
            group.sort_values("rank")["candidate_id"]
            .astype(str)
            .tolist()
        )
        if expected != actual:
            raise RuntimeError(
                f"Frozen candidate ranking mismatch for {arm}."
            )


def metric_delta(left, right):
    return {
        "spearman": float(
            left["spearman"] - right["spearman"]
        ),
        "top10_mean_percentile": float(
            left["top10_mean_percentile"]
            - right["top10_mean_percentile"]
        ),
        "top50_mean_percentile": float(
            left["top50_mean_percentile"]
            - right["top50_mean_percentile"]
        ),
        "top50_top1pct_hits": int(
            left["top50_top1pct_hits"]
            - right["top50_top1pct_hits"]
        ),
        "normalized_top1_regret": float(
            left["normalized_top1_regret"]
            - right["normalized_top1_regret"]
        ),
    }


def main(input_dir, stage_a_dir, out_dir):
    input_root = Path(input_dir)
    stage_a_root = Path(stage_a_dir)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    active_manifest = json.loads(
        (stage_a_root / "STAGE_A_MANIFEST.json").read_text(
            encoding="utf-8"
        )
    )

    predictions_path = (
        stage_a_root / "ACTIVE_BUDGET_PREDICTIONS.csv"
    )
    freeze_path = stage_a_root / "CANDIDATE_FREEZE.csv"
    acquisition_log_path = stage_a_root / "ACQUISITION_LOG.csv"
    hidden_truth_path = (
        input_root / ".HIDDEN_TRIPLE_TRUTH.sealed.csv"
    )

    checks = {
        "active_budget_predictions": (
            predictions_path,
            active_manifest["hashes"][
                "active_budget_predictions_sha256"
            ],
        ),
        "candidate_freeze": (
            freeze_path,
            active_manifest["hashes"][
                "candidate_freeze_sha256"
            ],
        ),
        "acquisition_log": (
            acquisition_log_path,
            active_manifest["hashes"][
                "acquisition_log_sha256"
            ],
        ),
        "sealed_hidden_truth": (
            hidden_truth_path,
            active_manifest["hashes"][
                "sealed_hidden_truth_sha256_expected"
            ],
        ),
    }

    verified = {}
    for name, (path, expected) in checks.items():
        actual = sha256_file(path)
        if actual != expected:
            raise RuntimeError(
                f"Hash mismatch for {name}: "
                f"{actual} != {expected}"
            )
        verified[name] = actual

    scores = pd.read_csv(predictions_path)
    freeze = pd.read_csv(freeze_path)
    verify_frozen_rankings(scores, freeze)

    passive_final = json.loads(
        PASSIVE_FINAL_PATH.read_text(encoding="utf-8")
    )

    # Hidden truth is loaded only after all active Stage-A hashes and
    # frozen candidate rankings verify.
    hidden_truth = pd.read_csv(hidden_truth_path)

    revealed = scores.merge(
        hidden_truth[["candidate_id", "DMS_score"]],
        on="candidate_id",
        how="inner",
        validate="one_to_one",
    )

    if len(revealed) != len(scores):
        raise RuntimeError(
            f"Truth coverage mismatch: "
            f"scores={len(scores)} revealed={len(revealed)}"
        )

    active_metrics = {}
    comparisons = {}
    threshold_vs_passive_40 = {}

    passive_40 = passive_final[
        "metrics_by_budget"
    ]["40"]["V8_3"]

    smallest_active_threshold = None

    summary_rows = []

    for budget in BUDGETS:
        prefix = f"A{budget:02d}"
        b3 = passive_eval.evaluate(
            revealed,
            f"{prefix}_B3",
        )
        b5 = passive_eval.evaluate(
            revealed,
            f"{prefix}_B5",
        )
        v83 = passive_eval.evaluate(
            revealed,
            f"{prefix}_V83",
        )

        active_metrics[str(budget)] = {
            "B3": b3,
            "B5": b5,
            "V8_3": v83,
        }

        passive_same = passive_final[
            "metrics_by_budget"
        ][str(budget)]["V8_3"]

        threshold = passive_eval.threshold_check(
            v83,
            passive_40,
        )

        threshold_vs_passive_40[str(budget)] = threshold

        if (
            smallest_active_threshold is None
            and threshold["passes_practical_threshold"]
        ):
            smallest_active_threshold = budget

        comparisons[str(budget)] = {
            "active_V8_3_vs_passive_same_budget": metric_delta(
                v83,
                passive_same,
            ),
            "active_V8_3_vs_passive_40pct": metric_delta(
                v83,
                passive_40,
            ),
        }

        summary_rows.append(
            {
                "budget_percent": budget,
                "visible_rows": active_manifest[
                    "diagnostics_by_budget"
                ][str(budget)]["visible_rows"],
                "router_mode": active_manifest[
                    "router_by_budget"
                ][str(budget)]["mode"],
                "active_spearman": v83["spearman"],
                "passive_same_budget_spearman": (
                    passive_same["spearman"]
                ),
                "delta_spearman_vs_passive_same": (
                    v83["spearman"]
                    - passive_same["spearman"]
                ),
                "active_top10_mean_percentile": (
                    v83["top10_mean_percentile"]
                ),
                "active_top50_mean_percentile": (
                    v83["top50_mean_percentile"]
                ),
                "active_top50_top1pct_hits": (
                    v83["top50_top1pct_hits"]
                ),
                "active_top1_regret": (
                    v83["normalized_top1_regret"]
                ),
                "passes_frozen_passive40_threshold": (
                    threshold["passes_practical_threshold"]
                ),
            }
        )

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(
        out / "ACTIVE_ACQUISITION_CURVE.csv",
        index=False,
    )

    if smallest_active_threshold is None:
        status = "ACTIVE_DID_NOT_REACH_PASSIVE40_THRESHOLD"
        reduction = None
        efficiency_factor = None
    else:
        status = "ACTIVE_REACHED_PASSIVE40_THRESHOLD"
        reduction = float(
            1.0 - (smallest_active_threshold / 40.0)
        )
        efficiency_factor = float(
            40.0 / smallest_active_threshold
        )

    final = {
        "version": "NABU_V8_3_RHLA_ACTIVE_ACQUISITION_FINAL_V1",
        "status": status,
        "development_only": True,
        "scientific_claim_boundary": (
            "RhlA truth was known from a prior experiment before this "
            "controller was designed. This result evaluates controller "
            "development retrospectively and is not fresh prospective proof."
        ),
        "hash_verification": verified,
        "fixed_hidden_rows": int(len(revealed)),
        "active_metrics_by_budget": active_metrics,
        "passive_reference_metrics_by_budget": {
            str(budget): passive_final[
                "metrics_by_budget"
            ][str(budget)]["V8_3"]
            for budget in BUDGETS
        },
        "comparisons": comparisons,
        "threshold_vs_frozen_passive_40pct": (
            threshold_vs_passive_40
        ),
        "smallest_active_budget_reaching_passive40_threshold_percent": (
            smallest_active_threshold
        ),
        "measurement_reduction_vs_passive_40pct": reduction,
        "measurement_efficiency_factor_vs_passive_40pct": (
            efficiency_factor
        ),
        "controller": {
            "exploration_weight": 0.5,
            "exploitation_weight": 0.5,
            "initial_budget_percent": 5,
            "active_targets_percent": [10, 20, 40],
        },
    }

    (out / "FINAL_RESULTS.json").write_text(
        json.dumps(final, indent=2),
        encoding="utf-8",
    )

    print(json.dumps(final, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        default="rhla_sample_efficiency_sealed_input",
    )
    parser.add_argument(
        "--stage-a",
        default="nabu_v8_3_rhla_active_acquisition_stage_a",
    )
    parser.add_argument(
        "--out",
        default="nabu_v8_3_rhla_active_acquisition_final",
    )
    args = parser.parse_args()
    main(args.input, args.stage_a, args.out)
