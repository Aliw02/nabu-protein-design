import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


BUDGETS = [5, 10, 20, 40]


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def percentile(values, x):
    values = np.asarray(values, dtype=float)
    return float(
        (np.sum(values < x) + 0.5 * np.sum(values == x))
        / len(values)
    )


def evaluate(frame, score_column):
    truth = frame["DMS_score"].to_numpy(dtype=float)
    ranked = frame.sort_values(
        [score_column, "candidate_id"],
        ascending=[False, True],
    )

    top1_true = float(ranked.iloc[0]["DMS_score"])
    top5 = ranked.head(5)["DMS_score"].to_numpy(dtype=float)
    top10 = ranked.head(10)["DMS_score"].to_numpy(dtype=float)
    top50 = ranked.head(50)["DMS_score"].to_numpy(dtype=float)

    cutoff1 = float(np.quantile(truth, 0.99))
    cutoff01 = float(np.quantile(truth, 0.999))
    best = float(np.max(truth))
    worst = float(np.min(truth))
    rho = float(
        spearmanr(
            frame[score_column].to_numpy(dtype=float),
            truth,
        ).statistic
    )

    return {
        "spearman": rho,
        "top1_true_score": top1_true,
        "top1_percentile": percentile(truth, top1_true),
        "top5_mean_percentile": float(
            np.mean([percentile(truth, x) for x in top5])
        ),
        "top10_mean_percentile": float(
            np.mean([percentile(truth, x) for x in top10])
        ),
        "top50_mean_percentile": float(
            np.mean([percentile(truth, x) for x in top50])
        ),
        "top10_top1pct_hits": int(np.sum(top10 >= cutoff1)),
        "top50_top1pct_hits": int(np.sum(top50 >= cutoff1)),
        "top10_top01pct_hits": int(np.sum(top10 >= cutoff01)),
        "top50_top01pct_hits": int(np.sum(top50 >= cutoff01)),
        "top50_mean_true_score": float(np.mean(top50)),
        "best_hidden_true_score": best,
        "normalized_top1_regret": (
            float((best - top1_true) / (best - worst))
            if best > worst
            else 0.0
        ),
    }


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


def safety_compare(router, b3):
    checks = {
        "spearman": router["spearman"] >= b3["spearman"],
        "top10_mean_percentile": (
            router["top10_mean_percentile"]
            >= b3["top10_mean_percentile"]
        ),
        "top50_mean_percentile": (
            router["top50_mean_percentile"]
            >= b3["top50_mean_percentile"]
        ),
        "top50_top1pct_hits": (
            router["top50_top1pct_hits"]
            >= b3["top50_top1pct_hits"]
        ),
        "normalized_top1_regret": (
            router["normalized_top1_regret"]
            <= b3["normalized_top1_regret"]
        ),
    }

    strict = {
        "spearman": router["spearman"] > b3["spearman"],
        "top10_mean_percentile": (
            router["top10_mean_percentile"]
            > b3["top10_mean_percentile"]
        ),
        "top50_mean_percentile": (
            router["top50_mean_percentile"]
            > b3["top50_mean_percentile"]
        ),
        "top50_top1pct_hits": (
            router["top50_top1pct_hits"]
            > b3["top50_top1pct_hits"]
        ),
        "normalized_top1_regret": (
            router["normalized_top1_regret"]
            < b3["normalized_top1_regret"]
        ),
    }

    return {
        "no_regression_by_metric": checks,
        "all_no_regression": bool(all(checks.values())),
        "strict_gain_by_metric": strict,
        "strict_gain_count": int(sum(strict.values())),
    }


def threshold_check(metric, ref):
    if ref["spearman"] > 0:
        spearman_retention = float(
            metric["spearman"] / ref["spearman"]
        )
    else:
        spearman_retention = None

    top10_gap = float(
        ref["top10_mean_percentile"]
        - metric["top10_mean_percentile"]
    )
    top50_gap = float(
        ref["top50_mean_percentile"]
        - metric["top50_mean_percentile"]
    )
    hits_loss = int(
        ref["top50_top1pct_hits"]
        - metric["top50_top1pct_hits"]
    )

    pass_spearman = (
        spearman_retention is not None
        and spearman_retention >= 0.90
    )
    pass_top10 = top10_gap <= 0.03
    pass_top50 = top50_gap <= 0.02
    pass_hits = hits_loss <= 1

    return {
        "spearman_retention_vs_40pct": spearman_retention,
        "top10_mean_percentile_gap_vs_40pct": top10_gap,
        "top50_mean_percentile_gap_vs_40pct": top50_gap,
        "top50_top1pct_hits_loss_vs_40pct": hits_loss,
        "pass_spearman_retention": bool(pass_spearman),
        "pass_top10_gap": bool(pass_top10),
        "pass_top50_gap": bool(pass_top50),
        "pass_hits_loss": bool(pass_hits),
        "passes_practical_threshold": bool(
            pass_spearman
            and pass_top10
            and pass_top50
            and pass_hits
        ),
    }


def main(input_dir, stage_a_dir, out_dir):
    input_root = Path(input_dir)
    stage_a_root = Path(stage_a_dir)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    manifest_path = stage_a_root / "STAGE_A_MANIFEST.json"
    manifest = json.loads(
        manifest_path.read_text(encoding="utf-8")
    )

    scores_path = stage_a_root / "ALL_BUDGET_PREDICTIONS.csv"
    freeze_path = stage_a_root / "CANDIDATE_FREEZE.csv"
    membership_path = stage_a_root / "BUDGET_MEMBERSHIP.csv"
    hidden_truth_path = input_root / ".HIDDEN_TRIPLE_TRUTH.sealed.csv"

    checks = {
        "budget_membership": (
            membership_path,
            manifest["hashes"]["budget_membership_sha256"],
        ),
        "all_budget_predictions": (
            scores_path,
            manifest["hashes"]["all_budget_predictions_sha256"],
        ),
        "candidate_freeze": (
            freeze_path,
            manifest["hashes"]["candidate_freeze_sha256"],
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

    scores = pd.read_csv(scores_path)
    freeze = pd.read_csv(freeze_path)
    verify_frozen_rankings(scores, freeze)

    # Hidden truth is loaded only after every Stage-A artifact verifies.
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

    metrics = {}
    safety = {}
    rows = []

    for percent in BUDGETS:
        prefix = f"P{percent:02d}"
        b3_col = f"{prefix}_B3"
        b5_col = f"{prefix}_B5"
        v83_col = f"{prefix}_V83"

        b3 = evaluate(revealed, b3_col)
        b5 = evaluate(revealed, b5_col)
        v83 = evaluate(revealed, v83_col)

        metrics[str(percent)] = {
            "B3": b3,
            "B5": b5,
            "V8_3": v83,
        }
        safety[str(percent)] = safety_compare(v83, b3)

        router_meta = manifest["router_by_budget"][str(percent)]
        diag = manifest["diagnostics_by_budget"][str(percent)]

        rows.append(
            {
                "budget_percent": percent,
                "visible_rows": diag["visible_rows"],
                "router_mode": router_meta["mode"],
                "B3_spearman": b3["spearman"],
                "V8_3_spearman": v83["spearman"],
                "delta_spearman": (
                    v83["spearman"] - b3["spearman"]
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
                "router_all_no_regression_vs_B3": (
                    safety[str(percent)]["all_no_regression"]
                ),
                "router_strict_gain_count_vs_B3": (
                    safety[str(percent)]["strict_gain_count"]
                ),
            }
        )

    controls = {
        "RANDOM_HASH": evaluate(revealed, "RANDOM_HASH"),
        "P40_SHUFFLE_V83": evaluate(
            revealed,
            "P40_SHUFFLE_V83",
        ),
    }

    reference = metrics["40"]["V8_3"]
    threshold = {}
    smallest = None

    for percent in BUDGETS:
        check = threshold_check(
            metrics[str(percent)]["V8_3"],
            reference,
        )
        threshold[str(percent)] = check
        if smallest is None and check["passes_practical_threshold"]:
            smallest = percent

    for row in rows:
        check = threshold[str(row["budget_percent"])]
        row.update(
            {
                "spearman_retention_vs_40pct": (
                    check["spearman_retention_vs_40pct"]
                ),
                "top10_gap_vs_40pct": (
                    check["top10_mean_percentile_gap_vs_40pct"]
                ),
                "top50_gap_vs_40pct": (
                    check["top50_mean_percentile_gap_vs_40pct"]
                ),
                "top50_hits_loss_vs_40pct": (
                    check["top50_top1pct_hits_loss_vs_40pct"]
                ),
                "passes_practical_threshold": (
                    check["passes_practical_threshold"]
                ),
            }
        )

    summary = pd.DataFrame(rows)
    summary.to_csv(
        out / "SAMPLE_EFFICIENCY_CURVE.csv",
        index=False,
    )

    final = {
        "version": "NABU_V8_3_RHLA_SAMPLE_EFFICIENCY_FINAL_V1",
        "status": (
            "PRACTICAL_THRESHOLD_FOUND"
            if smallest is not None
            else "NO_PRACTICAL_THRESHOLD_WITHIN_40PCT"
        ),
        "stage_A_manifest": manifest,
        "hash_verification": verified,
        "fixed_hidden_triple_rows": int(len(revealed)),
        "fixed_hidden_coverage_fraction": manifest[
            "fixed_hidden_coverage_fraction"
        ],
        "metrics_by_budget": metrics,
        "router_safety_by_budget": safety,
        "practical_threshold_checks": threshold,
        "smallest_practical_budget_percent": smallest,
        "reference_budget_percent": 40,
        "negative_controls": controls,
        "interpretation_boundary": (
            "All four nested budgets and all candidate rankings were "
            "frozen before hidden truth reveal. Results apply to the "
            "fixed hidden triple-mutant candidate pool scoreable by the 5% model."
        ),
    }

    (out / "FINAL_RESULTS.json").write_text(
        json.dumps(final, indent=2),
        encoding="utf-8",
    )

    revealed_rows = []
    columns = []
    for percent in BUDGETS:
        prefix = f"P{percent:02d}"
        columns.extend([
            (percent, f"{prefix}_B3"),
            (percent, f"{prefix}_B5"),
            (percent, f"{prefix}_V83"),
        ])
    columns.extend([
        (40, "RANDOM_HASH"),
        (40, "P40_SHUFFLE_V83"),
    ])

    for percent, score_column in columns:
        ranked = revealed.sort_values(
            [score_column, "candidate_id"],
            ascending=[False, True],
        ).head(50)

        for rank, row in enumerate(
            ranked.itertuples(),
            start=1,
        ):
            revealed_rows.append({
                "budget_percent": int(percent),
                "arm": score_column,
                "rank": int(rank),
                "candidate_id": row.candidate_id,
                "predicted_score": float(
                    getattr(row, score_column)
                ),
                "true_score": float(row.DMS_score),
            })

    pd.DataFrame(revealed_rows).to_csv(
        out / "REVEALED_TOP50.csv",
        index=False,
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
        default="nabu_v8_3_rhla_sample_efficiency_stage_a",
    )
    parser.add_argument(
        "--out",
        default="nabu_v8_3_rhla_sample_efficiency_final",
    )
    args = parser.parse_args()
    main(args.input, args.stage_a, args.out)
