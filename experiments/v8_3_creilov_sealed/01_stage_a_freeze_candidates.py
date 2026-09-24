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


def percentile(values, x):
    values = np.asarray(values, dtype=float)
    return float(
        (np.sum(values < x) + 0.5 * np.sum(values == x))
        / len(values)
    )


def rank_preserving_elite_rerank(frame, b3_col, b5_col, out_col, fraction=0.20):
    n = len(frame)
    k = max(1, int(math.ceil(n * fraction)))
    b3_order = frame.sort_values(
        [b3_col, "candidate_id"],
        ascending=[False, True],
    ).index.tolist()
    elite_idx = b3_order[:k]
    rest_idx = b3_order[k:]
    elite_sorted = (
        frame.loc[elite_idx]
        .sort_values(
            [b5_col, "candidate_id"],
            ascending=[False, True],
        )
        .index.tolist()
    )
    final_order = elite_sorted + rest_idx
    score = pd.Series(index=frame.index, dtype=float)
    for rank, idx in enumerate(final_order):
        score.loc[idx] = float(n - rank)
    frame[out_col] = score
    return {
        "elite_fraction": float(fraction),
        "elite_count": int(k),
        "boundary_preserved": True,
        "outside_B3_order_preserved": True,
    }


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

    elite_diag = oof_elite_diagnostic(
        visible_ids,
        visible_labels,
        model,
        fraction=0.20,
    )

    if global_delta > 0:
        frame[output_col] = frame[b5_col].astype(float)
        mode = "GLOBAL_HIGHER_ORDER"
        details = {"boundary_preserved": None}
    elif elite_diag["allow_elite_rerank"]:
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
            "reason": "Visible OOF did not support global or elite higher-order use.",
        }

    return {
        "B3_oof_spearman": b3_oof,
        "B4_oof_spearman": b4_oof,
        "B4_minus_B3": float(global_delta),
        "elite_diagnostic": elite_diag,
        "mode": mode,
        "details": details,
    }


def verify_seal(input_dir):
    root = Path(input_dir)
    manifest_path = root / "SEAL_MANIFEST.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    checks = {
        "visible": (
            root / "VISIBLE.csv",
            manifest["visible_sha256"],
        ),
        "hidden_ids": (
            root / "HIDDEN_IDS.csv",
            manifest["hidden_ids_sha256"],
        ),
    }
    for name, (path, expected) in checks.items():
        actual = sha256_file(path)
        if actual != expected:
            raise RuntimeError(
                f"Seal hash mismatch for {name}: {actual} != {expected}"
            )
    return manifest


def fit_and_score(visible, hidden, labels):
    visible_sets = visible["mutant"].map(v81.parse_mutations).tolist()
    visible_ids = visible["candidate_id"].astype(str).tolist()
    labels = np.asarray(labels, dtype=float)

    model = v81.fit_crossfitted_hierarchy(
        visible_sets,
        labels,
        visible_ids,
    )

    hidden = hidden.copy()
    hidden["mutation_set"] = hidden["mutant"].map(v81.parse_mutations)

    eligible = hidden[
        hidden["mutation_set"].map(
            lambda mutation_set: all(
                mutation in model["base"]["main"]
                for mutation in mutation_set
            )
            and any(
                pair in model["base"]["pair"]
                for pair in combinations(mutation_set, 2)
            )
        )
    ].copy()

    scored = [
        v81.score_hierarchy(mutation_set, model)
        for mutation_set in eligible["mutation_set"]
    ]
    scored_frame = pd.DataFrame(scored, index=eligible.index)
    for column in scored_frame.columns:
        eligible[column] = scored_frame[column]

    return model, eligible, visible_ids


def main(input_dir, out_dir):
    input_root = Path(input_dir)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    seal_manifest = verify_seal(input_root)

    visible = pd.read_csv(input_root / "VISIBLE.csv")
    hidden = pd.read_csv(input_root / "HIDDEN_IDS.csv")

    required_visible = {"candidate_id", "mutant", "DMS_score"}
    required_hidden = {"candidate_id", "mutant"}
    if not required_visible.issubset(visible.columns):
        raise RuntimeError(
            f"VISIBLE.csv missing columns: {required_visible - set(visible.columns)}"
        )
    if not required_hidden.issubset(hidden.columns):
        raise RuntimeError(
            f"HIDDEN_IDS.csv missing columns: {required_hidden - set(hidden.columns)}"
        )

    real_labels = visible["DMS_score"].to_numpy(dtype=float)
    real_model, real_scored, visible_ids = fit_and_score(
        visible,
        hidden,
        real_labels,
    )

    real_router = apply_router(
        real_scored,
        real_model,
        visible_ids,
        real_labels,
        "B3_RAW_PAIR",
        "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER",
        "V8_3_ADAPTIVE_ROUTER",
    )

    shuffled_labels = np.random.default_rng(SEED).permutation(real_labels)
    shuffle_model, shuffle_scored, shuffle_visible_ids = fit_and_score(
        visible,
        hidden,
        shuffled_labels,
    )

    common_ids = set(real_scored["candidate_id"]) & set(
        shuffle_scored["candidate_id"]
    )
    real_scored = real_scored[
        real_scored["candidate_id"].isin(common_ids)
    ].copy()
    shuffle_scored = shuffle_scored[
        shuffle_scored["candidate_id"].isin(common_ids)
    ].copy()

    shuffle_router = apply_router(
        shuffle_scored,
        shuffle_model,
        shuffle_visible_ids,
        shuffled_labels,
        "B3_RAW_PAIR",
        "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER",
        "SHUFFLED_LABEL_V8_3_CONTROL",
    )

    shuffle_scores = shuffle_scored[
        ["candidate_id", "SHUFFLED_LABEL_V8_3_CONTROL"]
    ].copy()

    scores = real_scored.merge(
        shuffle_scores,
        on="candidate_id",
        how="inner",
        validate="one_to_one",
    )

    keep_columns = [
        "candidate_id",
        "mutant",
        "mutation_count",
        "B3_RAW_PAIR",
        "B4_CROSSFIT_TRIPLET",
        "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER",
        "V8_3_ADAPTIVE_ROUTER",
        "SHUFFLED_LABEL_V8_3_CONTROL",
        "triplet_delta",
        "quartet_delta",
        "triplet_confidence",
        "quartet_confidence",
    ]
    keep_columns = [c for c in keep_columns if c in scores.columns]
    scores = scores[keep_columns].copy()

    scores_path = out / "ALL_HIDDEN_SCORES.csv"
    scores.to_csv(scores_path, index=False)

    arms = [
        "B3_RAW_PAIR",
        "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER",
        "V8_3_ADAPTIVE_ROUTER",
        "SHUFFLED_LABEL_V8_3_CONTROL",
    ]

    frozen_rows = []
    for arm in arms:
        ranked = scores.sort_values(
            [arm, "candidate_id"],
            ascending=[False, True],
        ).head(50)
        for rank, row in enumerate(ranked.itertuples(), start=1):
            frozen_rows.append(
                {
                    "arm": arm,
                    "rank": rank,
                    "candidate_id": row.candidate_id,
                    "mutation_count": int(row.mutation_count),
                    "predicted_score": float(getattr(row, arm)),
                    "budget10": bool(rank <= 10),
                    "budget50": True,
                }
            )

    freeze = pd.DataFrame(frozen_rows)
    freeze_path = out / "CANDIDATE_FREEZE.csv"
    freeze.to_csv(freeze_path, index=False)

    prereg_path = HERE.parent / "PREREGISTRATION.json"

    stage_a = {
        "version": "NABU_V8_3_CREILOV_STAGE_A_FREEZE_V1",
        "status": "CANDIDATES_FROZEN_BEFORE_HIDDEN_TRUTH_REVEAL",
        "base_freeze_commit": "c8afdcd7698231d95a39ef13a3fe22b6e3f507c5",
        "input_seal": seal_manifest,
        "visible_rows": int(len(visible)),
        "hidden_identity_rows": int(len(hidden)),
        "eligible_hidden_rows": int(len(scores)),
        "eligible_fraction": float(len(scores) / len(hidden)),
        "real_router": real_router,
        "shuffle_router": shuffle_router,
        "model_diagnostics": v81.model_diagnostics(real_model),
        "hashes": {
            "preregistration_sha256": sha256_file(prereg_path),
            "all_hidden_scores_sha256": sha256_file(scores_path),
            "candidate_freeze_sha256": sha256_file(freeze_path),
            "sealed_hidden_truth_sha256_expected": seal_manifest[
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
    print("STAGE A COMPLETE: candidate lists and hashes are frozen.")
    print("Do not modify Stage-A files before running Stage B.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        default="creilov_sealed_input",
    )
    parser.add_argument(
        "--out",
        default="nabu_v8_3_creilov_stage_a",
    )
    args = parser.parse_args()
    main(args.input, args.out)
