import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import pandas as pd


BUDGETS = [5, 10, 20, 40]
HERE = Path(__file__).resolve()
ROOT = HERE.parents[2]

EVAL_HELPER_PATH = (
    ROOT
    / "experiments"
    / "v8_3_rhla_sample_efficiency"
    / "02_stage_b_reveal_evaluate.py"
)

spec = importlib.util.spec_from_file_location(
    "nabu_eval_helper",
    EVAL_HELPER_PATH,
)
if spec is None or spec.loader is None:
    raise RuntimeError(
        f"Could not import frozen evaluator from {EVAL_HELPER_PATH}"
    )
eval_helper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(eval_helper)


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


def compare_primary(router, b3, tol=1e-12):
    no_regression = {
        "spearman": bool(
            router["spearman"] >= b3["spearman"] - tol
        ),
        "top10_mean_percentile": bool(
            router["top10_mean_percentile"]
            >= b3["top10_mean_percentile"] - tol
        ),
        "top50_mean_percentile": bool(
            router["top50_mean_percentile"]
            >= b3["top50_mean_percentile"] - tol
        ),
        "top50_top1pct_hits": bool(
            router["top50_top1pct_hits"]
            >= b3["top50_top1pct_hits"]
        ),
        "normalized_top1_regret": bool(
            router["normalized_top1_regret"]
            <= b3["normalized_top1_regret"] + tol
        ),
    }

    strict_gain = {
        "spearman": bool(
            router["spearman"] > b3["spearman"] + tol
        ),
        "top10_mean_percentile": bool(
            router["top10_mean_percentile"]
            > b3["top10_mean_percentile"] + tol
        ),
        "top50_mean_percentile": bool(
            router["top50_mean_percentile"]
            > b3["top50_mean_percentile"] + tol
        ),
        "top50_top1pct_hits": bool(
            router["top50_top1pct_hits"]
            > b3["top50_top1pct_hits"]
        ),
        "normalized_top1_regret": bool(
            router["normalized_top1_regret"]
            < b3["normalized_top1_regret"] - tol
        ),
    }

    return {
        "no_regression_by_metric": no_regression,
        "all_no_regression": bool(
            all(no_regression.values())
        ),
        "strict_gain_by_metric": strict_gain,
        "strict_gain_count": int(
            sum(strict_gain.values())
        ),
        "passes_core_gate": bool(
            all(no_regression.values())
            and sum(strict_gain.values()) >= 2
        ),
    }


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

    manifest = json.loads(
        (stage_a_root / "STAGE_A_MANIFEST.json").read_text(
            encoding="utf-8"
        )
    )

    predictions_path = (
        stage_a_root / "ALL_FINAL_GATE_PREDICTIONS.csv"
    )
    freeze_path = stage_a_root / "CANDIDATE_FREEZE.csv"
    passive_membership_path = (
        stage_a_root / "PASSIVE_BUDGET_MEMBERSHIP.csv"
    )
    acquisition_log_path = (
        stage_a_root / "ACTIVE_ACQUISITION_LOG.csv"
    )
    hidden_truth_path = (
        input_root / ".HIDDEN_4_5_TRUTH.sealed.csv"
    )

    checks = {
        "passive_membership": (
            passive_membership_path,
            manifest["hashes"][
                "passive_membership_sha256"
            ],
        ),
        "active_acquisition_log": (
            acquisition_log_path,
            manifest["hashes"][
                "active_acquisition_log_sha256"
            ],
        ),
        "all_predictions": (
            predictions_path,
            manifest["hashes"][
                "all_predictions_sha256"
            ],
        ),
        "candidate_freeze": (
            freeze_path,
            manifest["hashes"][
                "candidate_freeze_sha256"
            ],
        ),
        "sealed_hidden_truth": (
            hidden_truth_path,
            manifest["hashes"][
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

    verify_frozen_rankings(
        scores,
        freeze,
    )

    # Hidden truth is loaded exactly once, only after all Stage-A hashes
    # and candidate rankings have been verified.
    hidden_truth = pd.read_csv(
        hidden_truth_path
    )

    revealed = scores.merge(
        hidden_truth[
            ["candidate_id", "DMS_score"]
        ],
        on="candidate_id",
        how="inner",
        validate="one_to_one",
    )

    if len(revealed) != len(scores):
        raise RuntimeError(
            "Hidden truth coverage mismatch: "
            f"predictions={len(scores)} "
            f"revealed={len(revealed)}"
        )

    metrics = {
        "passive": {},
        "active": {},
    }

    summary_rows = []

    for family, prefix in [
        ("passive", "P"),
        ("active", "A"),
    ]:
        for budget in BUDGETS:
            arm_prefix = f"{prefix}{budget:02d}"

            b3 = eval_helper.evaluate(
                revealed,
                f"{arm_prefix}_B3",
            )
            b5 = eval_helper.evaluate(
                revealed,
                f"{arm_prefix}_B5",
            )
            v83 = eval_helper.evaluate(
                revealed,
                f"{arm_prefix}_V83",
            )

            metrics[family][str(budget)] = {
                "B3": b3,
                "B5": b5,
                "V8_3": v83,
            }

            summary_rows.append({
                "family": family,
                "budget_percent": int(budget),
                "router_mode": manifest[
                    "router_by_budget"
                ][family][str(budget)]["mode"],
                "B3_spearman": b3["spearman"],
                "V8_3_spearman": v83["spearman"],
                "delta_spearman_vs_B3": float(
                    v83["spearman"]
                    - b3["spearman"]
                ),
                "V8_3_top10_mean_percentile": (
                    v83["top10_mean_percentile"]
                ),
                "V8_3_top50_mean_percentile": (
                    v83["top50_mean_percentile"]
                ),
                "V8_3_top50_top1pct_hits": (
                    v83["top50_top1pct_hits"]
                ),
                "V8_3_top1_regret": (
                    v83["normalized_top1_regret"]
                ),
            })

    controls = {
        "RANDOM_HASH": eval_helper.evaluate(
            revealed,
            "RANDOM_HASH",
        ),
        "SHUFFLE40_V83": eval_helper.evaluate(
            revealed,
            "SHUFFLE40_V83",
        ),
    }

    core_gate = compare_primary(
        metrics["passive"]["40"]["V8_3"],
        metrics["passive"]["40"]["B3"],
    )

    passive_40_reference = metrics[
        "passive"
    ]["40"]["V8_3"]

    active_threshold = {}
    smallest_active_budget = None

    for budget in BUDGETS:
        threshold = eval_helper.threshold_check(
            metrics["active"][str(budget)]["V8_3"],
            passive_40_reference,
        )
        active_threshold[str(budget)] = threshold

        if (
            smallest_active_budget is None
            and threshold["passes_practical_threshold"]
        ):
            smallest_active_budget = budget

    active_target_pass = bool(
        active_threshold["20"][
            "passes_practical_threshold"
        ]
    )

    phase1_pass = bool(
        core_gate["passes_core_gate"]
        and active_target_pass
    )

    final_status = (
        "PHASE1_FINAL_GATE_PASS"
        if phase1_pass
        else "PHASE1_FINAL_GATE_CLOSED_WITH_LIMITATION"
    )

    comparisons = {
        "active_vs_passive_same_budget": {},
        "active_vs_passive40": {},
    }

    for budget in BUDGETS:
        comparisons[
            "active_vs_passive_same_budget"
        ][str(budget)] = metric_delta(
            metrics["active"][str(budget)][
                "V8_3"
            ],
            metrics["passive"][str(budget)][
                "V8_3"
            ],
        )

        comparisons[
            "active_vs_passive40"
        ][str(budget)] = metric_delta(
            metrics["active"][str(budget)][
                "V8_3"
            ],
            passive_40_reference,
        )

    summary = pd.DataFrame(
        summary_rows
    )
    summary.to_csv(
        out / "FINAL_GATE_CURVE.csv",
        index=False,
    )

    revealed_rows = []

    for family, prefix in [
        ("passive", "P"),
        ("active", "A"),
    ]:
        for budget in BUDGETS:
            for model_name in ["B3", "B5", "V83"]:
                arm = (
                    f"{prefix}{budget:02d}_"
                    f"{model_name}"
                )

                ranked = revealed.sort_values(
                    [arm, "candidate_id"],
                    ascending=[False, True],
                ).head(50)

                for rank, row in enumerate(
                    ranked.itertuples(),
                    start=1,
                ):
                    revealed_rows.append({
                        "family": family,
                        "budget_percent": int(budget),
                        "arm": model_name,
                        "rank": int(rank),
                        "candidate_id": row.candidate_id,
                        "mutation_count": int(
                            row.mutation_count
                        ),
                        "predicted_score": float(
                            getattr(row, arm)
                        ),
                        "true_score": float(
                            row.DMS_score
                        ),
                    })

    pd.DataFrame(
        revealed_rows
    ).to_csv(
        out / "REVEALED_TOP50.csv",
        index=False,
    )

    measurement_reduction = (
        None
        if smallest_active_budget is None
        else float(
            1.0
            - smallest_active_budget / 40.0
        )
    )

    efficiency_factor = (
        None
        if smallest_active_budget is None
        else float(
            40.0 / smallest_active_budget
        )
    )

    final = {
        "version": "NABU_V8_3_CR9114_H1_FINAL_PHASE1_GATE_RESULT_V1",
        "status": final_status,
        "phase1_closed_after_this_result": True,
        "no_post_reveal_retuning": True,
        "hash_verification": verified,
        "fixed_hidden_4_5_rows": int(
            len(revealed)
        ),
        "fixed_hidden_4_5_coverage_fraction": manifest[
            "fixed_common_hidden_coverage_fraction"
        ],
        "higher_order_qualification": manifest[
            "higher_order_qualification"
        ],
        "metrics": metrics,
        "negative_controls": controls,
        "core_v8_3_gate_passive40": core_gate,
        "active_threshold_vs_passive40": (
            active_threshold
        ),
        "active_20_target_pass": active_target_pass,
        "smallest_active_budget_reaching_passive40_threshold_percent": (
            smallest_active_budget
        ),
        "measurement_reduction_vs_passive40": (
            measurement_reduction
        ),
        "measurement_efficiency_factor_vs_passive40": (
            efficiency_factor
        ),
        "comparisons": comparisons,
        "final_gate": {
            "core_pass": bool(
                core_gate["passes_core_gate"]
            ),
            "active_20_pass": active_target_pass,
            "phase1_pass": phase1_pass,
        },
        "interpretation_boundary": (
            "One-shot fresh CR9114-H1 final Phase-1 gate on hidden "
            "4- and 5-mutation variants. Architecture and 50/50 active "
            "controller were frozen before this dataset outcome was revealed. "
            "Regardless of pass/fail, no further Phase-1 retuning is allowed; "
            "the next project step is Phase 2."
        ),
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
        default="cr9114_h1_final_gate_sealed_input",
    )
    parser.add_argument(
        "--stage-a",
        default="nabu_v8_3_cr9114_final_gate_stage_a",
    )
    parser.add_argument(
        "--out",
        default="nabu_v8_3_cr9114_final_gate_final",
    )
    args = parser.parse_args()
    main(args.input, args.stage_a, args.out)
