import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


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

    top1 = ranked.head(1)["DMS_score"].to_numpy(dtype=float)
    top5 = ranked.head(5)["DMS_score"].to_numpy(dtype=float)
    top10 = ranked.head(10)["DMS_score"].to_numpy(dtype=float)
    top50 = ranked.head(50)["DMS_score"].to_numpy(dtype=float)

    cutoff1 = float(np.quantile(truth, 0.99))
    cutoff01 = float(np.quantile(truth, 0.999))
    best = float(np.max(truth))
    worst = float(np.min(truth))
    top1_true = float(top1[0])

    rho = spearmanr(
        frame[score_column].to_numpy(dtype=float),
        truth,
    ).statistic

    return {
        "spearman": float(rho),
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


def replacement_diagnostic(frame, arm, baseline, k=50):
    arm_ids = set(
        frame.sort_values(
            [arm, "candidate_id"],
            ascending=[False, True],
        ).head(k)["candidate_id"]
    )
    base_ids = set(
        frame.sort_values(
            [baseline, "candidate_id"],
            ascending=[False, True],
        ).head(k)["candidate_id"]
    )

    added = arm_ids - base_ids
    removed = base_ids - arm_ids
    truth = frame.set_index("candidate_id")["DMS_score"].astype(float)

    added_values = (
        truth.loc[list(added)].to_numpy(dtype=float)
        if added
        else np.array([], dtype=float)
    )
    removed_values = (
        truth.loc[list(removed)].to_numpy(dtype=float)
        if removed
        else np.array([], dtype=float)
    )

    return {
        "overlap": int(len(arm_ids & base_ids)),
        "added_count": int(len(added)),
        "removed_count": int(len(removed)),
        "added_mean_true": (
            float(np.mean(added_values)) if len(added_values) else None
        ),
        "removed_mean_true": (
            float(np.mean(removed_values)) if len(removed_values) else None
        ),
        "replacement_gain_true": (
            float(np.mean(added_values) - np.mean(removed_values))
            if len(added_values) and len(removed_values)
            else None
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


def main(input_dir, stage_a_dir, out_dir):
    input_root = Path(input_dir)
    stage_a_root = Path(stage_a_dir)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    manifest_path = stage_a_root / "STAGE_A_MANIFEST.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    scores_path = stage_a_root / "ALL_HIDDEN_SCORES.csv"
    freeze_path = stage_a_root / "CANDIDATE_FREEZE.csv"
    hidden_truth_path = input_root / ".HIDDEN_TRUTH.sealed.csv"

    checks = {
        "all_hidden_scores": (
            scores_path,
            manifest["hashes"]["all_hidden_scores_sha256"],
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
    verified_hashes = {}
    for name, (path, expected) in checks.items():
        actual = sha256_file(path)
        if actual != expected:
            raise RuntimeError(
                f"Hash mismatch for {name}: {actual} != {expected}"
            )
        verified_hashes[name] = actual

    scores = pd.read_csv(scores_path)
    freeze = pd.read_csv(freeze_path)

    verify_frozen_rankings(scores, freeze)

    # Hidden truth is intentionally loaded only after every Stage-A hash
    # and frozen-ranking check has passed.
    hidden_truth = pd.read_csv(hidden_truth_path)

    revealed = scores.merge(
        hidden_truth[["candidate_id", "DMS_score"]],
        on="candidate_id",
        how="inner",
        validate="one_to_one",
    )

    if len(revealed) != len(scores):
        raise RuntimeError(
            f"Truth coverage mismatch: scores={len(scores)} "
            f"revealed={len(revealed)}"
        )

    arms = [
        "B3_RAW_PAIR",
        "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER",
        "V8_3_ADAPTIVE_ROUTER",
        "SHUFFLED_LABEL_V8_3_CONTROL",
    ]
    metrics = {
        arm: evaluate(revealed, arm)
        for arm in arms
    }

    b3 = metrics["B3_RAW_PAIR"]
    router = metrics["V8_3_ADAPTIVE_ROUTER"]
    shuffle = metrics["SHUFFLED_LABEL_V8_3_CONTROL"]

    tol = 1e-12
    primary = {
        "spearman": {
            "B3": b3["spearman"],
            "V8_3": router["spearman"],
            "no_regression": (
                router["spearman"] >= b3["spearman"] - tol
            ),
            "strict_gain": (
                router["spearman"] > b3["spearman"] + tol
            ),
        },
        "top10_mean_percentile": {
            "B3": b3["top10_mean_percentile"],
            "V8_3": router["top10_mean_percentile"],
            "no_regression": (
                router["top10_mean_percentile"]
                >= b3["top10_mean_percentile"] - tol
            ),
            "strict_gain": (
                router["top10_mean_percentile"]
                > b3["top10_mean_percentile"] + tol
            ),
        },
        "top50_mean_percentile": {
            "B3": b3["top50_mean_percentile"],
            "V8_3": router["top50_mean_percentile"],
            "no_regression": (
                router["top50_mean_percentile"]
                >= b3["top50_mean_percentile"] - tol
            ),
            "strict_gain": (
                router["top50_mean_percentile"]
                > b3["top50_mean_percentile"] + tol
            ),
        },
        "top50_top1pct_hits": {
            "B3": b3["top50_top1pct_hits"],
            "V8_3": router["top50_top1pct_hits"],
            "no_regression": (
                router["top50_top1pct_hits"]
                >= b3["top50_top1pct_hits"]
            ),
            "strict_gain": (
                router["top50_top1pct_hits"]
                > b3["top50_top1pct_hits"]
            ),
        },
        "normalized_top1_regret": {
            "B3": b3["normalized_top1_regret"],
            "V8_3": router["normalized_top1_regret"],
            "no_regression": (
                router["normalized_top1_regret"]
                <= b3["normalized_top1_regret"] + tol
            ),
            "strict_gain": (
                router["normalized_top1_regret"]
                < b3["normalized_top1_regret"] - tol
            ),
        },
    }

    all_no_regression = all(
        item["no_regression"] for item in primary.values()
    )
    strict_gain_count = sum(
        int(item["strict_gain"]) for item in primary.values()
    )
    promotion_pass = all_no_regression and strict_gain_count >= 2

    replacement = replacement_diagnostic(
        revealed,
        "V8_3_ADAPTIVE_ROUTER",
        "B3_RAW_PAIR",
        50,
    )

    result = {
        "version": "NABU_V8_3_EQFP611_SEALED_FINAL_V1",
        "status": (
            "SEALED_CANDIDATE_TEST_PASS"
            if promotion_pass
            else "SEALED_CANDIDATE_TEST_NO_PASS"
        ),
        "stage_A_manifest": manifest,
        "hash_verification": verified_hashes,
        "eligible_hidden_rows": int(len(revealed)),
        "metrics": metrics,
        "primary_gate": {
            "endpoints": primary,
            "all_no_regression": bool(all_no_regression),
            "strict_gain_count": int(strict_gain_count),
            "required_strict_gain_count": 2,
            "pass": bool(promotion_pass),
        },
        "V8_3_vs_B3": {
            "spearman": float(
                router["spearman"] - b3["spearman"]
            ),
            "top5_mean_percentile": float(
                router["top5_mean_percentile"]
                - b3["top5_mean_percentile"]
            ),
            "top10_mean_percentile": float(
                router["top10_mean_percentile"]
                - b3["top10_mean_percentile"]
            ),
            "top50_mean_percentile": float(
                router["top50_mean_percentile"]
                - b3["top50_mean_percentile"]
            ),
            "top50_top1pct_hits": int(
                router["top50_top1pct_hits"]
                - b3["top50_top1pct_hits"]
            ),
            "top50_top01pct_hits": int(
                router["top50_top01pct_hits"]
                - b3["top50_top01pct_hits"]
            ),
            "top50_mean_true_score": float(
                router["top50_mean_true_score"]
                - b3["top50_mean_true_score"]
            ),
            "normalized_top1_regret": float(
                router["normalized_top1_regret"]
                - b3["normalized_top1_regret"]
            ),
        },
        "top50_replacement_vs_B3": replacement,
        "negative_control": {
            "shuffle_spearman": shuffle["spearman"],
            "router_minus_shuffle_spearman": float(
                router["spearman"] - shuffle["spearman"]
            ),
            "router_top50_mean_percentile": router[
                "top50_mean_percentile"
            ],
            "shuffle_top50_mean_percentile": shuffle[
                "top50_mean_percentile"
            ],
        },
        "interpretation_boundary": (
            "One-shot sealed higher-order candidate-selection test on "
            "eqFP611 3-5 mutation variants with a deterministic identity-only "
            "70/30 split. Stage A required nonzero cross-fitted triplet and "
            "quartet memories before reveal. No architecture or threshold "
            "changes are permitted after this reveal."
        ),
    }

    result_path = out / "FINAL_RESULTS.json"
    result_path.write_text(
        json.dumps(result, indent=2),
        encoding="utf-8",
    )

    revealed_top = []
    for arm in arms:
        ranked = revealed.sort_values(
            [arm, "candidate_id"],
            ascending=[False, True],
        ).head(50)
        for rank, row in enumerate(ranked.itertuples(), start=1):
            revealed_top.append(
                {
                    "arm": arm,
                    "rank": rank,
                    "candidate_id": row.candidate_id,
                    "predicted_score": float(getattr(row, arm)),
                    "true_score": float(row.DMS_score),
                }
            )

    pd.DataFrame(revealed_top).to_csv(
        out / "REVEALED_TOP50.csv",
        index=False,
    )

    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        default="eqfp611_sealed_input",
    )
    parser.add_argument(
        "--stage-a",
        default="nabu_v8_3_eqfp611_stage_a",
    )
    parser.add_argument(
        "--out",
        default="nabu_v8_3_eqfp611_final",
    )
    args = parser.parse_args()
    main(args.input, args.stage_a, args.out)
