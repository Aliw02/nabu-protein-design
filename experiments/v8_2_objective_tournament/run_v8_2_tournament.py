import argparse
import importlib.util
import json
import math
import re
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
    raise RuntimeError(f"Could not import V8.1 runner from {V81_PATH}")
v81 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v81)


def safe_name(path):
    text = Path(path).stem
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text)


def pct(values, x):
    values = np.asarray(values, dtype=float)
    return float(
        (np.sum(values < x) + 0.5 * np.sum(values == x))
        / len(values)
    )


def evaluate(frame, score_column, truth_column):
    truth = frame[truth_column].to_numpy(dtype=float)
    ranked = frame.sort_values(
        [score_column, "candidate_id"],
        ascending=[False, True],
    )
    top1_cutoff = float(np.quantile(truth, 0.99))
    best = float(np.max(truth))
    worst = float(np.min(truth))

    def top_mean_percentile(k):
        values = ranked.head(min(k, len(ranked)))[truth_column].to_numpy(
            dtype=float
        )
        return float(np.mean([pct(truth, value) for value in values]))

    def top_hits(k):
        values = ranked.head(min(k, len(ranked)))[truth_column].to_numpy(
            dtype=float
        )
        return int(np.sum(values >= top1_cutoff))

    top1_truth = float(ranked.iloc[0][truth_column])
    top50 = ranked.head(min(50, len(ranked)))

    rho = spearmanr(
        frame[score_column].to_numpy(dtype=float),
        truth,
    ).statistic

    return {
        "spearman": float(rho),
        "top1_percentile": pct(truth, top1_truth),
        "top5_mean_percentile": top_mean_percentile(5),
        "top10_mean_percentile": top_mean_percentile(10),
        "top50_mean_percentile": top_mean_percentile(50),
        "top10_top1pct_hits": top_hits(10),
        "top50_top1pct_hits": top_hits(50),
        "top50_top1pct_enrichment": float(
            np.mean(
                top50[truth_column].to_numpy(dtype=float)
                >= top1_cutoff
            )
            / 0.01
        ),
        "top50_mean_true_score": float(
            top50[truth_column].astype(float).mean()
        ),
        "normalized_top1_regret": (
            float((best - top1_truth) / (best - worst))
            if best > worst
            else 0.0
        ),
    }


def top_ids(frame, score_column, k):
    return set(
        frame.sort_values(
            [score_column, "candidate_id"],
            ascending=[False, True],
        )
        .head(min(k, len(frame)))["candidate_id"]
        .tolist()
    )


def replacement_diagnostic(
    frame,
    score_column,
    baseline_column,
    truth_column,
    k,
):
    arm_ids = top_ids(frame, score_column, k)
    base_ids = top_ids(frame, baseline_column, k)
    added = arm_ids - base_ids
    displaced = base_ids - arm_ids

    truth_by_id = frame.set_index("candidate_id")[truth_column].astype(float)

    added_values = (
        truth_by_id.loc[list(added)].to_numpy(dtype=float)
        if added
        else np.array([], dtype=float)
    )
    displaced_values = (
        truth_by_id.loc[list(displaced)].to_numpy(dtype=float)
        if displaced
        else np.array([], dtype=float)
    )

    return {
        "overlap": int(len(arm_ids & base_ids)),
        "added_count": int(len(added)),
        "displaced_count": int(len(displaced)),
        "added_mean_true": (
            float(np.mean(added_values)) if len(added_values) else None
        ),
        "displaced_mean_true": (
            float(np.mean(displaced_values))
            if len(displaced_values)
            else None
        ),
        "replacement_gain_true": (
            float(np.mean(added_values) - np.mean(displaced_values))
            if len(added_values) and len(displaced_values)
            else None
        ),
    }


def rank_shift_diagnostic(frame, score_column, baseline_column):
    base_rank = (
        frame[baseline_column]
        .rank(method="average", ascending=False)
        .to_numpy(dtype=float)
    )
    arm_rank = (
        frame[score_column]
        .rank(method="average", ascending=False)
        .to_numpy(dtype=float)
    )
    shift = np.abs(arm_rank - base_rank)

    score_delta = (
        frame[score_column].to_numpy(dtype=float)
        - frame[baseline_column].to_numpy(dtype=float)
    )
    return {
        "mean_abs_rank_shift": float(np.mean(shift)),
        "p95_abs_rank_shift": float(np.quantile(shift, 0.95)),
        "mean_abs_score_delta": float(np.mean(np.abs(score_delta))),
        "p95_abs_score_delta": float(
            np.quantile(np.abs(score_delta), 0.95)
        ),
    }


def add_variant_scores(frame, oof_prefers_higher_order):
    frame = frame.copy()
    b3 = frame["B3_RAW_PAIR"].to_numpy(dtype=float)
    b5 = frame[
        "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER"
    ].to_numpy(dtype=float)
    triplet_delta = frame["triplet_delta"].to_numpy(dtype=float)
    quartet_delta = frame["quartet_delta"].to_numpy(dtype=float)
    triplet_conf = frame["triplet_confidence"].to_numpy(dtype=float)
    quartet_conf = frame["quartet_confidence"].to_numpy(dtype=float)

    full_delta = b5 - b3

    # Direction A: higher-order dosage.
    frame["T_HALF_HIGHER_ORDER"] = b3 + 0.50 * full_delta
    frame["T_QUARTER_HIGHER_ORDER"] = b3 + 0.25 * full_delta
    frame["T_THREE_QUARTER_HIGHER_ORDER"] = b3 + 0.75 * full_delta

    # Direction B: a second confidence gate. V8.1 deltas are already
    # confidence-aware; this deliberately tests whether sparse landscapes
    # need stronger shrinkage.
    frame["T_CONFIDENCE_SQUARED"] = (
        b3
        + triplet_delta * triplet_conf
        + quartet_delta * quartet_conf
    )

    # Direction C: objective split. Let higher-order reorder only the
    # B3 elite region, preserving most of the whole-landscape ranking.
    for fraction in (0.10, 0.20, 0.30):
        cutoff = float(np.quantile(b3, 1.0 - fraction))
        mask = b3 >= cutoff
        name = f"T_ELITE_{int(fraction * 100)}PCT_HIGHER_ORDER"
        score = b3.copy()
        score[mask] = b5[mask]
        frame[name] = score

    # Direction D: scale-free rank fusion between pairwise and V8.1.
    b3_rank = frame["B3_RAW_PAIR"].rank(
        method="average", pct=True, ascending=True
    )
    b5_rank = frame[
        "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER"
    ].rank(method="average", pct=True, ascending=True)
    frame["T_RANK_FUSION_50_50"] = (
        0.50 * b3_rank + 0.50 * b5_rank
    )
    frame["T_RANK_FUSION_PAIR_75"] = (
        0.75 * b3_rank + 0.25 * b5_rank
    )

    # Direction E: visible-OOF landscape gate. This never reads hidden truth
    # to decide whether HO is global or elite-only.
    if oof_prefers_higher_order:
        frame["T_OOF_OBJECTIVE_GATE"] = b5
    else:
        cutoff = float(np.quantile(b3, 0.80))
        mask = b3 >= cutoff
        score = b3.copy()
        score[mask] = b5[mask]
        frame["T_OOF_OBJECTIVE_GATE"] = score

    return frame


def prepare_dataset(csv_path):
    data = pd.read_csv(csv_path)
    mutation_column = v81.find_column(
        data.columns,
        v81.MUTATION_COLUMNS,
    )
    truth_column = v81.find_column(
        data.columns,
        v81.FITNESS_COLUMNS,
    )

    if mutation_column is None or truth_column is None:
        raise ValueError(
            f"{csv_path}: need mutation and fitness columns; "
            f"found {list(data.columns)}"
        )

    data = data[[mutation_column, truth_column]].dropna().copy()
    data["mutation_set"] = data[mutation_column].map(v81.parse_mutations)
    data = data[
        data["mutation_set"].map(
            lambda x: 3 <= len(x) <= 5 and "*" not in "".join(x)
        )
    ].copy()
    data["candidate_id"] = data["mutation_set"].map(
        lambda x: ":".join(x)
    )
    data["bucket"] = data["candidate_id"].map(
        lambda x: v81.fnv1a32(x) % 10
    )
    data["mutation_count"] = data["mutation_set"].map(len)

    visible = data[data["bucket"] <= 6].copy()
    hidden = data[data["bucket"] >= 7].copy()

    visible_sets = visible["mutation_set"].tolist()
    visible_labels = visible[truth_column].to_numpy(dtype=float)
    visible_ids = visible["candidate_id"].tolist()

    model = v81.fit_crossfitted_hierarchy(
        visible_sets,
        visible_labels,
        visible_ids,
    )

    eligible = hidden[
        hidden["mutation_set"].map(
            lambda mutation_set: all(
                mutation in model["base"]["main"]
                for mutation in mutation_set
            )
            and any(
                pair in model["base"]["pair"]
                for pair in __import__("itertools").combinations(
                    mutation_set, 2
                )
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

    oof_b3_rho = float(
        spearmanr(model["oof_b3"], visible_labels).statistic
    )
    oof_b4_rho = float(
        spearmanr(model["oof_b4"], visible_labels).statistic
    )
    oof_delta = oof_b4_rho - oof_b3_rho

    eligible = add_variant_scores(
        eligible,
        oof_prefers_higher_order=oof_delta > 0.0,
    )

    return {
        "data": data,
        "visible": visible,
        "hidden": hidden,
        "eligible": eligible,
        "truth_column": truth_column,
        "model": model,
        "oof_b3_rho": oof_b3_rho,
        "oof_b4_rho": oof_b4_rho,
        "oof_delta": float(oof_delta),
    }


def analyze_dataset(csv_path, out_root):
    prepared = prepare_dataset(csv_path)
    eligible = prepared["eligible"]
    truth_column = prepared["truth_column"]
    dataset_name = safe_name(csv_path)

    dataset_out = out_root / dataset_name
    dataset_out.mkdir(parents=True, exist_ok=True)

    arms = [
        "B3_RAW_PAIR",
        "B4_CROSSFIT_TRIPLET",
        "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER",
        "T_QUARTER_HIGHER_ORDER",
        "T_HALF_HIGHER_ORDER",
        "T_THREE_QUARTER_HIGHER_ORDER",
        "T_CONFIDENCE_SQUARED",
        "T_ELITE_10PCT_HIGHER_ORDER",
        "T_ELITE_20PCT_HIGHER_ORDER",
        "T_ELITE_30PCT_HIGHER_ORDER",
        "T_RANK_FUSION_50_50",
        "T_RANK_FUSION_PAIR_75",
        "T_OOF_OBJECTIVE_GATE",
    ]

    metrics = {
        arm: evaluate(eligible, arm, truth_column)
        for arm in arms
    }
    baseline = metrics["B3_RAW_PAIR"]

    per_arm = {}
    rows = []

    for arm in arms:
        m = metrics[arm]
        replacements = {
            str(k): replacement_diagnostic(
                eligible,
                arm,
                "B3_RAW_PAIR",
                truth_column,
                k,
            )
            for k in (5, 10, 50)
        }
        shifts = rank_shift_diagnostic(
            eligible,
            arm,
            "B3_RAW_PAIR",
        )

        deltas = {
            "spearman": float(
                m["spearman"] - baseline["spearman"]
            ),
            "top5_mean_percentile": float(
                m["top5_mean_percentile"]
                - baseline["top5_mean_percentile"]
            ),
            "top10_mean_percentile": float(
                m["top10_mean_percentile"]
                - baseline["top10_mean_percentile"]
            ),
            "top50_mean_percentile": float(
                m["top50_mean_percentile"]
                - baseline["top50_mean_percentile"]
            ),
            "top50_top1pct_hits": int(
                m["top50_top1pct_hits"]
                - baseline["top50_top1pct_hits"]
            ),
            "top50_mean_true_score": float(
                m["top50_mean_true_score"]
                - baseline["top50_mean_true_score"]
            ),
            "normalized_top1_regret": float(
                m["normalized_top1_regret"]
                - baseline["normalized_top1_regret"]
            ),
        }

        tags = []
        if deltas["spearman"] > 0:
            tags.append("GLOBAL_RANK_GAIN")
        elif deltas["spearman"] < 0:
            tags.append("GLOBAL_RANK_HARM")

        elite_gain = (
            deltas["top5_mean_percentile"] > 0
            or deltas["top50_mean_percentile"] > 0
            or deltas["top50_top1pct_hits"] > 0
        )
        tags.append("ELITE_GAIN" if elite_gain else "NO_ELITE_GAIN")

        rep50 = replacements["50"]["replacement_gain_true"]
        if rep50 is not None:
            tags.append(
                "TOP50_REPLACEMENTS_BETTER"
                if rep50 > 0
                else "TOP50_REPLACEMENTS_WORSE"
            )

        per_arm[arm] = {
            "metrics": m,
            "deltas_vs_B3": deltas,
            "replacement_diagnostics": replacements,
            "rank_shift_diagnostic": shifts,
            "diagnostic_tags": tags,
        }

        rows.append(
            {
                "dataset": dataset_name,
                "arm": arm,
                **m,
                "delta_spearman_vs_B3": deltas["spearman"],
                "delta_top5_pct_vs_B3": deltas[
                    "top5_mean_percentile"
                ],
                "delta_top50_pct_vs_B3": deltas[
                    "top50_mean_percentile"
                ],
                "delta_top50_hits_vs_B3": deltas[
                    "top50_top1pct_hits"
                ],
                "delta_top50_true_vs_B3": deltas[
                    "top50_mean_true_score"
                ],
                "top50_overlap_with_B3": replacements["50"][
                    "overlap"
                ],
                "top50_replacement_gain_true": replacements["50"][
                    "replacement_gain_true"
                ],
                "mean_abs_rank_shift": shifts[
                    "mean_abs_rank_shift"
                ],
            }
        )

    result = {
        "version": "NABU_V8_2_DEVELOPMENT_TOURNAMENT",
        "status": "REVEALED_DEVELOPMENT_ONLY_NOT_PROSPECTIVE",
        "dataset": Path(csv_path).name,
        "rows": int(len(prepared["data"])),
        "visible": int(len(prepared["visible"])),
        "hidden": int(len(prepared["hidden"])),
        "eligible_hidden": int(len(eligible)),
        "visible_oof_landscape_diagnostic": {
            "B3_spearman": prepared["oof_b3_rho"],
            "B4_spearman": prepared["oof_b4_rho"],
            "B4_minus_B3": prepared["oof_delta"],
            "OOF_gate_mode": (
                "GLOBAL_HIGHER_ORDER"
                if prepared["oof_delta"] > 0
                else "ELITE_20PCT_HIGHER_ORDER"
            ),
        },
        "arms": per_arm,
        "model_diagnostics": v81.model_diagnostics(
            prepared["model"]
        ),
        "interpretation_boundary": (
            "This tournament intentionally uses already revealed development "
            "landscapes to improve NABU. It is for architecture diagnosis, "
            "not a fresh blind validation claim."
        ),
    }

    (dataset_out / "V8_2_DATASET_RESULTS.json").write_text(
        json.dumps(result, indent=2),
        encoding="utf-8",
    )
    eligible.drop(columns=["mutation_set"]).to_csv(
        dataset_out / "V8_2_CANDIDATE_SCORES.csv",
        index=False,
    )
    pd.DataFrame(rows).to_csv(
        dataset_out / "V8_2_ARM_SUMMARY.csv",
        index=False,
    )

    return result, rows


def aggregate(all_results, all_rows):
    frame = pd.DataFrame(all_rows)
    if frame.empty:
        return {}, frame

    baseline_name = "B3_RAW_PAIR"
    arms = [
        arm for arm in frame["arm"].unique().tolist()
        if arm != baseline_name
    ]

    aggregate = {}
    for arm in arms:
        group = frame[frame["arm"] == arm].copy()
        aggregate[arm] = {
            "datasets": int(len(group)),
            "mean_delta_spearman_vs_B3": float(
                group["delta_spearman_vs_B3"].mean()
            ),
            "mean_delta_top5_pct_vs_B3": float(
                group["delta_top5_pct_vs_B3"].mean()
            ),
            "mean_delta_top50_pct_vs_B3": float(
                group["delta_top50_pct_vs_B3"].mean()
            ),
            "sum_delta_top50_hits_vs_B3": int(
                group["delta_top50_hits_vs_B3"].sum()
            ),
            "mean_delta_top50_true_vs_B3": float(
                group["delta_top50_true_vs_B3"].mean()
            ),
            "global_rank_gain_datasets": int(
                np.sum(group["delta_spearman_vs_B3"] > 0)
            ),
            "top5_gain_datasets": int(
                np.sum(group["delta_top5_pct_vs_B3"] > 0)
            ),
            "top50_gain_datasets": int(
                np.sum(group["delta_top50_pct_vs_B3"] > 0)
            ),
            "positive_top50_replacement_datasets": int(
                np.sum(
                    group["top50_replacement_gain_true"].fillna(0) > 0
                )
            ),
        }

    return aggregate, frame


def main(csvs, out_dir):
    out_root = Path(out_dir)
    out_root.mkdir(parents=True, exist_ok=True)

    all_results = {}
    all_rows = []

    for csv_path in csvs:
        print(f"[V8.2] dataset={csv_path}")
        result, rows = analyze_dataset(csv_path, out_root)
        all_results[result["dataset"]] = result
        all_rows.extend(rows)

    aggregate_summary, aggregate_frame = aggregate(
        all_results,
        all_rows,
    )

    final = {
        "version": "NABU_V8_2_MULTI_LANDSCAPE_DEVELOPMENT_TOURNAMENT",
        "status": "REVEALED_DEVELOPMENT_ONLY_NOT_PROSPECTIVE",
        "datasets": list(all_results.keys()),
        "aggregate_vs_B3": aggregate_summary,
        "selection_rule": (
            "Do not choose from one scalar score. Prefer arms that retain "
            "or improve elite selection across landscapes while limiting "
            "whole-landscape rank harm; inspect replacement diagnostics "
            "to understand why each arm succeeds or fails."
        ),
    }

    (out_root / "V8_2_TOURNAMENT_SUMMARY.json").write_text(
        json.dumps(final, indent=2),
        encoding="utf-8",
    )
    aggregate_frame.to_csv(
        out_root / "V8_2_TOURNAMENT_ALL_ARMS.csv",
        index=False,
    )

    print(json.dumps(final, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "csvs",
        nargs="+",
        help="One or more V8.1-compatible landscape CSV files.",
    )
    parser.add_argument(
        "--out",
        default="nabu_v8_2_tournament_results",
    )
    args = parser.parse_args()
    main(args.csvs, args.out)
