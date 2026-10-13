"""Open external targets only after validating the complete prediction freeze."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from common import HERE, read_json, sha256, source_manifest, write_json
from metrics import metrics, paired_bootstrap


def evaluate(out):
    if (out / "DECISION.json").exists():
        raise RuntimeError("Evaluation already exists; inspect preserved results instead of overwriting")
    freeze = read_json(out / "PREDICTION_FREEZE.json")
    assert source_manifest() == freeze["source_sha256"], "Code changed after prediction freeze"
    for name, expected in freeze["artifacts_sha256"].items():
        assert sha256(out / name) == expected, name
    replay_check = read_json(out / "DETERMINISM.json")
    assert replay_check["passed"] and replay_check["tuning_choices_and_scores_identical"]
    assert replay_check["original_prediction_freeze_sha256"] == sha256(out / "PREDICTION_FREEZE.json")
    for name, expected in replay_check["identical_artifact_sha256"].items():
        assert sha256(out / name) == expected
    manifest = read_json(out / "PRE_RUN_MANIFEST.json")
    assert sha256(HERE / "DATA_AUDIT.json") == manifest["data_audit_sha256"]
    for name, expected in manifest["prepared_sha256"].items():
        assert sha256(HERE / "data/prepared" / name) == expected, name
    config = read_json(HERE / "protocol.json")
    methods = [config["primary_candidate"], *config["baselines"], *config["controls"]]
    selection = read_json(out / "SELECTION.json")
    champion, candidate = selection["baseline_champion"], config["primary_candidate"]
    assert champion == freeze["baseline_champion"]
    assert freeze["external_targets_used"] is False
    dev = pd.read_csv(HERE / "data/prepared/development.csv")
    ext_ids = pd.read_csv(HERE / "data/prepared/external_identities.csv")
    dev_pred = pd.read_csv(out / "development_predictions.csv")
    ext_pred = pd.read_csv(out / "external_predictions.csv")
    for data, prediction in [(dev, dev_pred), (ext_ids, ext_pred)]:
        assert set(prediction.method) == set(methods)
        assert len(prediction) == len(data) * len(methods)
        assert not prediction.duplicated(["candidate_id", "method"]).any()
        assert np.isfinite(prediction.prediction_pEC50).all()
        for method in methods:
            part = prediction[prediction.method == method]
            assert set(zip(part.candidate_id, part.sequence)) == set(zip(data.candidate_id, data.sequence))
    # This is the first evaluation-stage read of the external target values.
    truth = pd.read_csv(HERE / "data/prepared/external_truth.csv")
    assert len(truth) == len(ext_ids)
    dev_cases = dev.merge(dev_pred.drop(columns="outer_fold"), on=["candidate_id", "sequence"], validate="one_to_many")
    ext_cases = truth.merge(ext_pred.drop(columns="outer_fold"), on=["candidate_id", "sequence"], validate="one_to_many")
    assert len(dev_cases) == 125 * len(methods) and len(ext_cases) == 19 * len(methods)
    equivalent = ext_cases.groupby(["method", "sequence"], sort=True).agg(
        target_pEC50=("target_pEC50", "mean"), prediction_pEC50=("prediction_pEC50", "mean"),
        candidate_id=("candidate_id", lambda values: "|".join(sorted(values))), record_count=("candidate_id", "size"))
    equivalent = equivalent.reset_index()
    assert equivalent.groupby("method").size().eq(18).all()
    duplicate_prediction_spread = ext_cases.groupby(["method", "sequence"]).prediction_pEC50.agg(lambda x: x.max() - x.min())
    assert duplicate_prediction_spread.max() < 1e-12
    all_metrics, rows = {}, []
    partitions = {"development_oof": dev_cases, "external_sequence_units": equivalent, "external_all_records": ext_cases}
    for partition, data in partitions.items():
        all_metrics[partition] = {}
        for method in methods:
            frame = data[data.method == method].sort_values("candidate_id")
            result = metrics(frame.target_pEC50, frame.prediction_pEC50, frame.candidate_id.tolist())
            active = frame[frame.target_pEC50 > 9.0]
            result["better_than_1nM_subset"] = (metrics(active.target_pEC50, active.prediction_pEC50, active.candidate_id.tolist())
                                              if len(active) >= 5 else {"n": len(active), "metrics": None})
            all_metrics[partition][method] = result
            rows.append({"partition": partition, "method": method, **{k: v for k, v in result.items() if k not in ("top5_ids", "better_than_1nM_subset")}})
    # Recompute the champion and dev scores from the frozen per-case artifact.
    champion_check = min(config["baselines"], key=lambda m: (all_metrics["development_oof"][m]["rmse"], m))
    assert champion == champion_check
    for method in methods:
        assert abs(all_metrics["development_oof"][method]["rmse"] - selection["development_metrics"][method]["rmse"]) < 1e-10
    fold_metrics = []
    for (fold, method), data in dev_cases.groupby(["outer_fold", "method"]):
        fold_metrics.append({"outer_fold": int(fold), "method": method,
                             **metrics(data.target_pEC50, data.prediction_pEC50, data.candidate_id.tolist())})
    bootstraps = {}
    for index, partition in enumerate(["development_oof", "external_sequence_units"]):
        data = partitions[partition]
        left = data[data.method == candidate].sort_values("candidate_id")
        right = data[data.method == champion].set_index("candidate_id").loc[left.candidate_id]
        group = left.sequence_group.to_numpy() if partition == "development_oof" else np.arange(len(left))
        bootstraps[partition] = paired_bootstrap(left.target_pEC50.to_numpy(), left.prediction_pEC50.to_numpy(),
                                               right.prediction_pEC50.to_numpy(), group,
                                               config["bootstrap_resamples"], config["seed"] + index)
    d, e = all_metrics["development_oof"], all_metrics["external_sequence_units"]
    factor = 1 - config["gate"]["minimum_rmse_relative_gain"]
    def rank_at_least(a, b, gain=0.0):
        return a is not None and b is not None and a >= b + gain
    checks = {
        "integrity_and_full_coverage": True,
        "unchanged_deterministic_replay": True,
        "development_rmse_gain_at_least_5pct": d[candidate]["rmse"] <= factor * d[champion]["rmse"],
        "development_spearman_not_lower": rank_at_least(d[candidate]["spearman"], d[champion]["spearman"]),
        "development_beats_additive": d[candidate]["rmse"] < d["JOINT_ADDITIVE"]["rmse"],
        "development_beats_shuffle": d[candidate]["rmse"] < d["SHUFFLED_JOINT_PAIRS"]["rmse"],
        "external_rmse_gain_at_least_5pct": e[candidate]["rmse"] <= factor * e[champion]["rmse"],
        "external_spearman_gain_at_least_005": rank_at_least(e[candidate]["spearman"], e[champion]["spearman"],
                                                            config["gate"]["minimum_external_spearman_gain"]),
        "external_top5_mean_not_lower": e[candidate]["top5_mean_measured_pEC50"] >= e[champion]["top5_mean_measured_pEC50"],
        "external_beats_additive": e[candidate]["rmse"] < e["JOINT_ADDITIVE"]["rmse"],
        "development_paired_mse_ci_upper_negative": bootstraps["development_oof"]["upper_bound_negative"],
        "external_paired_mse_ci_upper_negative": bootstraps["external_sequence_units"]["upper_bound_negative"],
    }
    decision = {"experiment": config["experiment"],
                "decision": "ADVANCE_TO_INDEPENDENT_REPLICATION" if all(checks.values()) else "DO_NOT_ADVANCE",
                "evaluated_utc": datetime.now(timezone.utc).isoformat(),
                "baseline_champion_chosen_on_development": champion, "candidate": candidate,
                "checks": {k: bool(v) for k, v in checks.items()}, "failed_checks": [k for k, v in checks.items() if not v],
                "prediction_freeze_sha256": sha256(out / "PREDICTION_FREEZE.json"),
                "development_relative_rmse_gain": 1 - d[candidate]["rmse"] / d[champion]["rmse"],
                "external_relative_rmse_gain": 1 - e[candidate]["rmse"] / e[champion]["rmse"],
                "claim_scope": "Retrospective supervised GLP1R pEC50 prediction only; no de novo or wet-lab claim"}
    duplicate_truth = truth[truth.sequence.duplicated(False)].copy()
    write_json(out / "external_duplicate_sensitivity.json", {
        "records": duplicate_truth.to_dict("records"),
        "target_pEC50_range": float(duplicate_truth.target_pEC50.max() - duplicate_truth.target_pEC50.min()),
        "interpretation": "Sequence-equivalent records; chemical identity/replicate status not established",
        "identical_sequence_prediction_max_spread": float(duplicate_prediction_spread.max())})
    for name, frame in [("development_per_case.csv", dev_cases), ("external_per_record.csv", ext_cases),
                        ("external_per_sequence.csv", equivalent)]:
        frame["error_pEC50"] = frame.prediction_pEC50 - frame.target_pEC50
        frame["squared_error"] = frame.error_pEC50 ** 2
        frame.to_csv(out / name, index=False)
    pd.DataFrame(rows).to_csv(out / "metrics.csv", index=False)
    write_json(out / "metrics.json", all_metrics)
    write_json(out / "fold_metrics.json", fold_metrics)
    write_json(out / "paired_bootstrap.json", bootstraps)
    write_json(out / "DECISION.json", decision)
    write_report(out, all_metrics, bootstraps, decision)
    print(pd.DataFrame(rows)[["partition", "method", "rmse", "spearman", "top5_mean_measured_pEC50"]].to_string(index=False))
    print(decision["decision"], decision["failed_checks"])


def write_report(out, results, bootstraps, decision):
    candidate, champion = decision["candidate"], decision["baseline_champion_chosen_on_development"]
    lines = ["# GLP1R activity bridge: frozen result", "", f"**Decision: {decision['decision']}**", "",
             f"Development-selected baseline: `{champion}`. Candidate: `{candidate}`.", "",
             "All models used the same lossless 44-slot sequence representation. Development used five outer edit-distance-grouped folds with three inner grouped folds; external predictions were frozen before evaluation. External primary metrics use 18 unique sequences from 19 published records.", "",
             "| Method | Development RMSE | Development Spearman | External RMSE | External Spearman | External top-5 mean pEC50 |",
             "|---|---:|---:|---:|---:|---:|"]
    def number(value):
        return "undefined (constant)" if value is None else f"{value:.6f}"
    for method, dev in results["development_oof"].items():
        ext = results["external_sequence_units"][method]
        lines.append(f"| {method} | {number(dev['rmse'])} | {number(dev['spearman'])} | {number(ext['rmse'])} | {number(ext['spearman'])} | {number(ext['top5_mean_measured_pEC50'])} |")
    lines += ["", "Lower RMSE is better; higher Spearman and top-5 mean potency are better. These are predictions of a published cell-assay potency endpoint, not therapeutic efficacy, molecular binding validation, or new peptide experiments.", "", "## Frozen gate", ""]
    lines += [f"- {'PASS' if passed else 'FAIL'}: {name}" for name, passed in decision["checks"].items()]
    lines += ["", "## Paired uncertainty", "", "Candidate-minus-champion MSE; negative values favor the candidate.", ""]
    for partition, values in bootstraps.items():
        lines.append(f"- {partition}: {values['candidate_minus_baseline_mse']:.6f}; 95% CI [{values['ci95_low']:.6f}, {values['ci95_high']:.6f}], {values['resampling_units']} resampling units.")
    lines += ["", "## Boundaries", "",
              "This is a standard categorical ANOVA ridge model testing a NABU restructuring hypothesis, not a novel algorithm. The legacy arms apply the original B2/B3 equations to the new adapter and are not a rerun of the complete V8.3 hierarchy. All training activity labels are allowed here, changing the question from the historical low-order-only transfer task.", "",
              "Public retrospective data, one small related peptide family, incomplete chemical metadata, study-to-study assay differences, and sequence-equivalent external records limit generalization. Grouped development validation reduces near-clone leakage; it does not establish protein-family transfer. The published paper uses other protocols, so these results are not a like-for-like victory over that paper.", "",
              "No post-outcome retuning is permitted. The per-case artifacts, raw source hashes, complete grids, support diagnostics, freeze records, all 19-record sensitivity metrics, and negative controls are retained. NucB and the canonical model are untouched. Competition eligibility and de novo functional validation remain separate."]
    (out / "REPORT.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=HERE / "results")
    args = parser.parse_args()
    try:
        evaluate(args.output.resolve())
    except Exception as error:
        # Never overwrite an existing scientific decision, even on an attempted rerun.
        if not (args.output / "DECISION.json").exists():
            write_json(args.output / "DECISION.json", {"decision": "INVALID", "error": repr(error)})
        raise
