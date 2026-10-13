"""Fit the preregistered arms and freeze predictions before external evaluation."""
from __future__ import annotations

import os
for variable in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[variable] = "1"

import argparse
from datetime import datetime, timezone
import platform
import subprocess
from pathlib import Path
from time import perf_counter

import numpy as np
import pandas as pd
import scipy
import sklearn
from threadpoolctl import threadpool_info, threadpool_limits

from common import HERE, codes_from, read_json, sha256, source_manifest, write_json
from metrics import metrics
from models import fit_predict


def git(*arguments):
    return subprocess.check_output(["git", *arguments], cwd=HERE, text=True).strip()


def predict(out):
    config = read_json(HERE / "protocol.json")
    audit = read_json(HERE / "DATA_AUDIT.json")
    if out.exists() and any(out.iterdir()):
        raise RuntimeError("Refusing to overwrite an existing run; use a new --output for exact replay")
    for name, expected in config["source_sha256"].items():
        assert sha256(HERE / "data/raw" / name) == expected
    for name, expected in audit["prepared_sha256"].items():
        assert sha256(HERE / "data/prepared" / name) == expected
    sources = source_manifest()
    frozen_paths = [*sources, "DATA_AUDIT.json", *["data/prepared/" + name for name in audit["prepared_sha256"]],
                    *["data/raw/" + name for name in config["source_sha256"]]]
    for path in frozen_paths:
        git("ls-files", "--error-unmatch", "--", path)
    if git("diff", "HEAD", "--", *frozen_paths):
        raise RuntimeError("Experiment sources and inputs must be committed before fitting")
    assert read_json(HERE / "INTEGRITY_TESTS.json")["passed"] is True
    started = perf_counter()
    out.mkdir(parents=True)
    write_json(out / "PRE_RUN_MANIFEST.json", {
        "started_utc": datetime.now(timezone.utc).isoformat(), "git_commit": git("rev-parse", "HEAD"),
        "source_sha256": sources, "prepared_sha256": audit["prepared_sha256"],
        "protocol_sha256": sha256(HERE / "protocol.json"), "data_audit_sha256": sha256(HERE / "DATA_AUDIT.json"),
        "integrity_tests_sha256": sha256(HERE / "INTEGRITY_TESTS.json"),
        "external_targets_read_by_predictor": False})
    write_json(out / "software_versions.json", {
        "python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
        "scipy": scipy.__version__, "scikit_learn": sklearn.__version__, "platform": platform.platform(),
        "threadpools": threadpool_info()})
    dev = pd.read_csv(HERE / "data/prepared/development.csv")
    external = pd.read_csv(HERE / "data/prepared/external_identities.csv")
    assert "target_pEC50" not in external
    methods = [config["primary_candidate"], *config["baselines"], *config["controls"]]
    x, query = codes_from(dev), codes_from(external)
    length_columns = ["prefix_length", "core_length", "suffix_length"]
    lengths, query_lengths = dev[length_columns].to_numpy(), external[length_columns].to_numpy()
    y, groups = dev.target_pEC50.to_numpy(), dev.sequence_group.to_numpy()
    dev_rows, external_rows, supports, tuning = [], [], [], []
    fold_ids = dev.outer_fold.to_numpy()

    def prediction_rows(indices, source, predictions, method, partition, fold, train_x, query_x, support):
        records = []
        known = [set(train_x[:, j]) for j in range(train_x.shape[1])]
        for position, index in enumerate(indices):
            record = {"candidate_id": source.iloc[index].candidate_id,
                      "sequence": source.iloc[index].sequence, "method": method,
                      "prediction_pEC50": float(predictions[position]), "outer_fold": fold}
            records.append(record)
            if method == config["primary_candidate"] or support is not None:
                diagnostics = {"partition": partition, "outer_fold": fold, "method": method,
                               "candidate_id": record["candidate_id"],
                               "unseen_position_states": sum(state not in known[j] for j, state in enumerate(query_x[position])),
                               "feature_width": train_x.shape[1]}
                if support is not None:
                    diagnostics.update(support[position])
                supports.append(diagnostics)
        return records

    with threadpool_limits(limits=1):
        for fold in sorted(np.unique(fold_ids)):
            fit, hold = np.flatnonzero(fold_ids != fold), np.flatnonzero(fold_ids == fold)
            assert not set(groups[fit]).intersection(groups[hold])
            for method in methods:
                predictions, trace, support = fit_predict(
                    method, x[fit], y[fit], lengths[fit], groups[fit], x[hold], lengths[hold],
                    config, config["seed"] + int(fold))
                tuning.append({"partition": "development", "outer_fold": int(fold), "method": method, **trace})
                dev_rows.extend(prediction_rows(hold, dev, predictions, method, "development", int(fold),
                                                x[fit], x[hold], support))
            print(f"Completed grouped development fold {fold + 1}/{config['outer_folds']}", flush=True)

        # The champion is selected before any external label is loaded.
        development = pd.DataFrame(dev_rows)
        dev_metrics = {}
        for method in methods:
            aligned = dev.merge(development[development.method == method], on=["candidate_id", "sequence"], validate="one_to_one")
            dev_metrics[method] = metrics(aligned.target_pEC50, aligned.prediction_pEC50, aligned.candidate_id.tolist())
        champion = min(config["baselines"], key=lambda name: (dev_metrics[name]["rmse"], name))
        write_json(out / "SELECTION.json", {"baseline_champion": champion, "selection_metric": "development_pooled_oof_rmse",
                                            "development_metrics": dev_metrics, "external_targets_used": False})
        for method in methods:
            predictions, trace, support = fit_predict(method, x, y, lengths, groups, query, query_lengths,
                                                     config, config["seed"] + 1000)
            tuning.append({"partition": "external", "outer_fold": -1, "method": method, **trace})
            external_rows.extend(prediction_rows(range(len(external)), external, predictions, method,
                                                 "external", -1, x, query, support))
        print("Completed external predictions without loading external targets", flush=True)
    development.to_csv(out / "development_predictions.csv", index=False)
    pd.DataFrame(external_rows).to_csv(out / "external_predictions.csv", index=False)
    pd.DataFrame(supports).to_csv(out / "support_diagnostics.csv", index=False)
    write_json(out / "tuning_trace.json", tuning)
    write_json(out / "PREDICTION_FREEZE.json", {
        "frozen_utc": datetime.now(timezone.utc).isoformat(), "git_commit": git("rev-parse", "HEAD"),
        "prediction_stage_seconds": perf_counter() - started,
        "baseline_champion": champion, "external_targets_used": False,
        "artifacts_sha256": {p.name: sha256(p) for p in sorted(out.iterdir()) if p.is_file()},
        "source_sha256": sources})
    print(f"Prediction freeze saved: {out / 'PREDICTION_FREEZE.json'}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=HERE / "results")
    predict(parser.parse_args().output.resolve())
