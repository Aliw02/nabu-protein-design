import argparse
import hashlib
import importlib.util
import json
import math
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


BUDGETS = [5, 10, 20, 40]
SEED = 161

HERE = Path(__file__).resolve()
V81_PATH = (
    HERE.parents[1]
    / "v8_1_crossfit_higher_order_dev"
    / "run_v8_1_crossfit_dev.py"
)

spec = importlib.util.spec_from_file_location("nabu_v81", V81_PATH)
if spec is None or spec.loader is None:
    raise RuntimeError(f"Could not import frozen V8.1 runner from {V81_PATH}")
v81 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v81)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def fnv1a32(text):
    h = 2166136261
    for byte in str(text).encode("utf-8"):
        h ^= byte
        h = (h * 16777619) & 0xFFFFFFFF
    return h


def percentile(values, x):
    values = np.asarray(values, dtype=float)
    return float(
        (np.sum(values < x) + 0.5 * np.sum(values == x))
        / len(values)
    )


def oof_elite_diagnostic(candidate_ids, labels, model, fraction=0.20):
    work = pd.DataFrame(
        {
            "candidate_id": list(candidate_ids),
            "truth": np.asarray(labels, dtype=float),
            "B3_OOF": np.asarray(model["oof_b3"], dtype=float),
            "B4_OOF": np.asarray(model["oof_b4"], dtype=float),
        }
    )

    n = len(work)
    k_region = max(1, int(math.ceil(n * fraction)))
    region = (
        work.sort_values(
            ["B3_OOF", "candidate_id"],
            ascending=[False, True],
        )
        .head(k_region)
        .copy()
    )

    top_k = min(50, len(region))
    b3_top = region.sort_values(
        ["B3_OOF", "candidate_id"],
        ascending=[False, True],
    ).head(top_k)
    b4_top = region.sort_values(
        ["B4_OOF", "candidate_id"],
        ascending=[False, True],
    ).head(top_k)

    all_truth = work["truth"].to_numpy(dtype=float)
    cutoff1 = float(np.quantile(all_truth, 0.99))

    def summarize(rows):
        vals = rows["truth"].to_numpy(dtype=float)
        return {
            "top50_mean_percentile": float(
                np.mean([percentile(all_truth, x) for x in vals])
            ),
            "top50_top1pct_hits": int(np.sum(vals >= cutoff1)),
            "top50_mean_visible_score": float(np.mean(vals)),
        }

    b3 = summarize(b3_top)
    b4 = summarize(b4_top)
    allow = (
        b4["top50_mean_percentile"] > b3["top50_mean_percentile"]
        and b4["top50_top1pct_hits"] >= b3["top50_top1pct_hits"]
    )

    return {
        "region_fraction": float(fraction),
        "region_count": int(k_region),
        "B3": b3,
        "B4": b4,
        "delta_top50_mean_percentile": float(
            b4["top50_mean_percentile"]
            - b3["top50_mean_percentile"]
        ),
        "delta_top50_top1pct_hits": int(
            b4["top50_top1pct_hits"]
            - b3["top50_top1pct_hits"]
        ),
        "allow_elite_rerank": bool(allow),
    }


def rank_preserving_elite_rerank(
    frame,
    b3_col,
    b5_col,
    output_col,
    fraction=0.20,
):
    n = len(frame)
    k = max(1, int(math.ceil(n * fraction)))

    b3_order = frame.sort_values(
        [b3_col, "candidate_id"],
        ascending=[False, True],
    ).index.tolist()
    elite = b3_order[:k]
    rest = b3_order[k:]

    elite_sorted = (
        frame.loc[elite]
        .sort_values(
            [b5_col, "candidate_id"],
            ascending=[False, True],
        )
        .index.tolist()
    )

    order = elite_sorted + rest
    score = pd.Series(index=frame.index, dtype=float)
    for rank, idx in enumerate(order):
        score.loc[idx] = float(n - rank)
    frame[output_col] = score

    return {
        "elite_fraction": float(fraction),
        "elite_count": int(k),
        "boundary_preserved": True,
        "outside_B3_order_preserved": True,
    }


def apply_router(
    frame,
    model,
    visible_ids,
    visible_labels,
    b3_col,
    b5_col,
    output_col,
):
    b3_oof = float(
        spearmanr(model["oof_b3"], visible_labels).statistic
    )
    b4_oof = float(
        spearmanr(model["oof_b4"], visible_labels).statistic
    )
    global_delta = b4_oof - b3_oof

    elite = oof_elite_diagnostic(
        visible_ids,
        visible_labels,
        model,
        fraction=0.20,
    )

    if global_delta > 0:
        frame[output_col] = frame[b5_col].astype(float)
        mode = "GLOBAL_HIGHER_ORDER"
        details = {"boundary_preserved": None}
    elif elite["allow_elite_rerank"]:
        details = rank_preserving_elite_rerank(
            frame,
            b3_col,
            b5_col,
            output_col,
            fraction=0.20,
        )
        mode = "RANK_PRESERVING_B3_TOP20_B5_RERANK"
    else:
        frame[output_col] = frame[b3_col].astype(float)
        mode = "B3_PROTECTED_NO_HIGHER_ORDER"
        details = {
            "boundary_preserved": True,
            "reason": (
                "Visible OOF did not support global or elite "
                "higher-order use."
            ),
        }

    return {
        "B3_oof_spearman": b3_oof,
        "B4_oof_spearman": b4_oof,
        "B4_minus_B3": float(global_delta),
        "elite_diagnostic": elite,
        "mode": mode,
        "details": details,
    }


def verify_seal(input_dir):
    root = Path(input_dir)
    manifest = json.loads(
        (root / "SEAL_MANIFEST.json").read_text(encoding="utf-8")
    )

    checks = {
        "training_pool": (
            root / "TRAINING_POOL.csv",
            manifest["training_pool_sha256"],
        ),
        "hidden_ids": (
            root / "HIDDEN_TRIPLE_IDS.csv",
            manifest["hidden_ids_sha256"],
        ),
    }

    for name, (path, expected) in checks.items():
        actual = sha256_file(path)
        if actual != expected:
            raise RuntimeError(
                f"Seal hash mismatch for {name}: "
                f"{actual} != {expected}"
            )

    return manifest


def fit_model(visible):
    mutation_sets = visible["mutant"].map(v81.parse_mutations).tolist()
    labels = visible["DMS_score"].to_numpy(dtype=float)
    ids = visible["candidate_id"].astype(str).tolist()

    model = v81.fit_crossfitted_hierarchy(
        mutation_sets,
        labels,
        ids,
    )
    return model, ids, labels


def eligible_mask(hidden, model):
    mutation_sets = hidden["mutant"].map(v81.parse_mutations)
    return mutation_sets.map(
        lambda mutation_set: all(
            mutation in model["base"]["main"]
            for mutation in mutation_set
        )
        and any(
            pair in model["base"]["pair"]
            for pair in combinations(mutation_set, 2)
        )
    )


def score_candidates(hidden, model):
    scored = hidden.copy()
    scored["mutation_set"] = scored["mutant"].map(v81.parse_mutations)

    rows = [
        v81.score_hierarchy(mutation_set, model)
        for mutation_set in scored["mutation_set"]
    ]
    sf = pd.DataFrame(rows, index=scored.index)

    for column in sf.columns:
        scored[column] = sf[column]

    return scored


def budget_count(total_eligible, percent):
    return int(math.ceil(total_eligible * (percent / 100.0)))


def main(input_dir, out_dir):
    input_root = Path(input_dir)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    seal = verify_seal(input_root)
    training_pool = pd.read_csv(input_root / "TRAINING_POOL.csv")
    hidden_ids = pd.read_csv(input_root / "HIDDEN_TRIPLE_IDS.csv")

    required_train = {
        "candidate_id",
        "mutant",
        "DMS_score",
        "mutation_count",
    }
    required_hidden = {
        "candidate_id",
        "mutant",
        "mutation_count",
    }
    if not required_train.issubset(training_pool.columns):
        raise RuntimeError(
            "TRAINING_POOL.csv missing columns: "
            f"{required_train - set(training_pool.columns)}"
        )
    if not required_hidden.issubset(hidden_ids.columns):
        raise RuntimeError(
            "HIDDEN_IDS.csv missing columns: "
            f"{required_hidden - set(hidden_ids.columns)}"
        )

    total_eligible = int(seal["total_eligible_rows"])

    ordered = training_pool.copy()
    ordered["budget_hash"] = ordered["candidate_id"].map(
        lambda cid: fnv1a32("BUDGET|" + str(cid))
    )
    ordered = ordered.sort_values(
        ["budget_hash", "candidate_id"],
        ascending=[True, True],
    ).reset_index(drop=True)

    counts = {
        percent: budget_count(total_eligible, percent)
        for percent in BUDGETS
    }
    if max(counts.values()) > len(ordered):
        raise RuntimeError(
            "Training pool is smaller than the preregistered 40% "
            f"budget: training_pool={len(ordered)} "
            f"required={max(counts.values())}"
        )

    budget_frames = {
        percent: ordered.head(counts[percent]).copy()
        for percent in BUDGETS
    }

    membership_rows = []
    for row in ordered.itertuples():
        record = {
            "candidate_id": row.candidate_id,
            "budget_hash": int(row.budget_hash),
        }
        for percent in BUDGETS:
            record[f"in_{percent:02d}pct"] = bool(
                row.Index < counts[percent]
            )
        membership_rows.append(record)

    membership = pd.DataFrame(membership_rows)
    membership_path = out / "BUDGET_MEMBERSHIP.csv"
    membership.to_csv(membership_path, index=False)

    first_budget = BUDGETS[0]
    model_5, ids_5, labels_5 = fit_model(
        budget_frames[first_budget]
    )
    mask_5 = eligible_mask(hidden_ids, model_5)
    fixed_hidden = hidden_ids[mask_5].copy()

    coverage = float(len(fixed_hidden) / len(hidden_ids))
    if coverage < 0.50 or len(fixed_hidden) < 50:
        abort = {
            "status": "ABORT_BEFORE_REVEAL",
            "reason": (
                "5% model does not cover enough of the fixed hidden "
                "candidate pool."
            ),
            "hidden_rows": int(len(hidden_ids)),
            "fixed_hidden_rows": int(len(fixed_hidden)),
            "coverage_fraction": coverage,
            "required_coverage_fraction": 0.50,
        }
        (out / "STAGE_A_QUALIFICATION_ABORT.json").write_text(
            json.dumps(abort, indent=2),
            encoding="utf-8",
        )
        print(json.dumps(abort, indent=2))
        raise SystemExit(2)

    base_scores = fixed_hidden[
        ["candidate_id", "mutant", "mutation_count"]
    ].copy()

    diagnostics_by_budget = {}
    router_by_budget = {}

    for percent in BUDGETS:
        visible = budget_frames[percent]
        model, visible_ids, labels = fit_model(visible)

        # Because budgets are nested, every candidate scoreable at 5%
        # must remain scoreable later. Abort if that invariant fails.
        current_mask = eligible_mask(fixed_hidden, model)
        if not bool(current_mask.all()):
            raise RuntimeError(
                f"Nested scoreability invariant failed at {percent}%."
            )

        scored = score_candidates(fixed_hidden, model)

        prefix = f"P{percent:02d}"
        b3_col = f"{prefix}_B3"
        b4_col = f"{prefix}_B4"
        b5_col = f"{prefix}_B5"
        v83_col = f"{prefix}_V83"

        scored[b3_col] = scored["B3_RAW_PAIR"].astype(float)
        scored[b4_col] = scored["B4_CROSSFIT_TRIPLET"].astype(float)
        scored[b5_col] = scored[
            "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER"
        ].astype(float)

        router = apply_router(
            scored,
            model,
            visible_ids,
            labels,
            b3_col,
            b5_col,
            v83_col,
        )

        keep = scored[
            ["candidate_id", b3_col, b4_col, b5_col, v83_col]
        ].copy()

        base_scores = base_scores.merge(
            keep,
            on="candidate_id",
            how="inner",
            validate="one_to_one",
        )

        diag = v81.model_diagnostics(model)
        diagnostics_by_budget[str(percent)] = {
            "visible_rows": int(len(visible)),
            "visible_fraction_of_all_eligible": float(
                len(visible) / total_eligible
            ),
            "model_diagnostics": diag,
        }
        router_by_budget[str(percent)] = router

    base_scores["RANDOM_HASH"] = base_scores["candidate_id"].map(
        lambda cid: float(
            fnv1a32("RANDOM|" + str(cid)) / 4294967296.0
        )
    )

    visible_40 = budget_frames[40]
    labels_40 = visible_40["DMS_score"].to_numpy(dtype=float)
    shuffled = np.random.default_rng(SEED).permutation(labels_40)

    visible_40_shuffle = visible_40.copy()
    visible_40_shuffle["DMS_score"] = shuffled
    shuffle_model, shuffle_ids, shuffle_labels = fit_model(
        visible_40_shuffle
    )
    shuffle_scored = score_candidates(fixed_hidden, shuffle_model)
    shuffle_scored["P40_SHUFFLE_B3"] = shuffle_scored[
        "B3_RAW_PAIR"
    ].astype(float)
    shuffle_scored["P40_SHUFFLE_B5"] = shuffle_scored[
        "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER"
    ].astype(float)

    shuffle_router = apply_router(
        shuffle_scored,
        shuffle_model,
        shuffle_ids,
        shuffle_labels,
        "P40_SHUFFLE_B3",
        "P40_SHUFFLE_B5",
        "P40_SHUFFLE_V83",
    )

    base_scores = base_scores.merge(
        shuffle_scored[
            ["candidate_id", "P40_SHUFFLE_V83"]
        ],
        on="candidate_id",
        how="inner",
        validate="one_to_one",
    )

    scores_path = out / "ALL_BUDGET_PREDICTIONS.csv"
    base_scores.to_csv(scores_path, index=False)

    freeze_rows = []
    score_columns = []
    for percent in BUDGETS:
        prefix = f"P{percent:02d}"
        score_columns.extend(
            [
                f"{prefix}_B3",
                f"{prefix}_B5",
                f"{prefix}_V83",
            ]
        )
    score_columns.extend(["RANDOM_HASH", "P40_SHUFFLE_V83"])

    for score_column in score_columns:
        ranked = base_scores.sort_values(
            [score_column, "candidate_id"],
            ascending=[False, True],
        ).head(50)

        for rank, row in enumerate(ranked.itertuples(), start=1):
            freeze_rows.append(
                {
                    "arm": score_column,
                    "rank": rank,
                    "candidate_id": row.candidate_id,
                    "mutation_count": int(row.mutation_count),
                    "predicted_score": float(
                        getattr(row, score_column)
                    ),
                    "budget10": bool(rank <= 10),
                    "budget50": True,
                }
            )

    freeze = pd.DataFrame(freeze_rows)
    freeze_path = out / "CANDIDATE_FREEZE.csv"
    freeze.to_csv(freeze_path, index=False)

    prereg_path = HERE.parent / "PREREGISTRATION.json"

    stage_a = {
        "version": "NABU_V8_3_RHLA_SAMPLE_EFFICIENCY_STAGE_A_V1",
        "status": "ALL_BUDGETS_FROZEN_BEFORE_HIDDEN_TRUTH_REVEAL",
        "base_freeze_commit": (
            "c8afdcd7698231d95a39ef13a3fe22b6e3f507c5"
        ),
        "input_seal": seal,
        "total_eligible_rows": int(total_eligible),
        "training_pool_rows": int(len(training_pool)),
        "hidden_triple_identity_rows": int(len(hidden_ids)),
        "fixed_hidden_rows": int(len(fixed_hidden)),
        "fixed_hidden_coverage_fraction": coverage,
        "budget_counts": {
            str(k): int(v) for k, v in counts.items()
        },
        "router_by_budget": router_by_budget,
        "diagnostics_by_budget": diagnostics_by_budget,
        "negative_controls": {
            "random_hash": "deterministic FNV1a32 ranking",
            "P40_SHUFFLE_V83_router": shuffle_router,
            "P40_SHUFFLE_V83_model_diagnostics": (
                v81.model_diagnostics(shuffle_model)
            ),
        },
        "hashes": {
            "preregistration_sha256": sha256_file(prereg_path),
            "budget_membership_sha256": sha256_file(membership_path),
            "all_budget_predictions_sha256": sha256_file(scores_path),
            "candidate_freeze_sha256": sha256_file(freeze_path),
            "sealed_hidden_truth_sha256_expected": seal[
                "hidden_truth_sha256"
            ],
        },
        "hidden_truth_loaded": False,
    }

    manifest_path = out / "STAGE_A_MANIFEST.json"
    manifest_path.write_text(
        json.dumps(stage_a, indent=2),
        encoding="utf-8",
    )

    print(json.dumps(stage_a, indent=2))
    print("STAGE A COMPLETE: all budgets frozen before reveal.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        default="rhla_sample_efficiency_sealed_input",
    )
    parser.add_argument(
        "--out",
        default="nabu_v8_3_rhla_sample_efficiency_stage_a",
    )
    args = parser.parse_args()
    main(args.input, args.out)
