"""Validate unchanged CPU replay before any external outcome evaluation."""
from __future__ import annotations

import argparse
from pathlib import Path

from common import HERE, read_json, sha256, write_json


def check(original, replay):
    first = read_json(original / "PREDICTION_FREEZE.json")
    second = read_json(replay / "PREDICTION_FREEZE.json")
    assert first["source_sha256"] == second["source_sha256"]
    names = ["development_predictions.csv", "external_predictions.csv", "SELECTION.json", "support_diagnostics.csv"]
    hashes = {}
    for name in names:
        hashes[name] = sha256(original / name)
        assert hashes[name] == sha256(replay / name), name
    def without_timings(path):
        return [{k: v for k, v in entry.items() if k != "fit_and_predict_seconds"}
                for entry in read_json(path / "tuning_trace.json")]
    assert without_timings(original) == without_timings(replay)
    path = original / "DETERMINISM.json"
    if path.exists():
        raise RuntimeError("Replay verification already recorded")
    write_json(path, {"passed": True, "identical_artifact_sha256": hashes,
                      "tuning_choices_and_scores_identical": True,
                      "original_prediction_freeze_sha256": sha256(original / "PREDICTION_FREEZE.json"),
                      "replay_prediction_freeze_sha256": sha256(replay / "PREDICTION_FREEZE.json"),
                      "replay_prediction_seconds": second["prediction_stage_seconds"],
                      "external_labels_used_for_check": False})
    print("PASS: predictions, support, champion selection and all inner tuning scores are identical")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--original", type=Path, default=HERE / "results")
    parser.add_argument("--replay", type=Path, default=HERE / "replay")
    args = parser.parse_args()
    check(args.original, args.replay)
