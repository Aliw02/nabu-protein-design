from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

EXPECTED_MODEL_REVISION = "c731040fcd8d73dceaa04b0a8e6329b345b0f5df"
EXPECTED_MODEL_SHA256 = (
    "24c5fa474c48f3b754b86efe752d5f189d2bcd88190fa2270fc92b2ef3034189"
)

ARMS = ["LOW_ONLY", "REPRESENTATION_ONLY", "FULL", "SHUFFLED_FULL"]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def finite_spearman(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) < 2 or np.std(a) == 0 or np.std(b) == 0:
        return None
    x = float(spearmanr(a, b).statistic)
    return x if np.isfinite(x) else None


def finite_pearson(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) < 2 or np.std(a) == 0 or np.std(b) == 0:
        return None
    x = float(pearsonr(a, b).statistic)
    return x if np.isfinite(x) else None


def regression_metrics(target, prediction):
    target = np.asarray(target, dtype=float)
    prediction = np.asarray(prediction, dtype=float)
    rmse = float(np.sqrt(np.mean((target - prediction) ** 2)))
    return {
        "count": int(len(target)),
        "spearman": finite_spearman(target, prediction),
        "pearson": finite_pearson(target, prediction),
        "rmse": rmse,
        "target_std": float(np.std(target)),
    }


def rank_metrics(target, prediction):
    return {
        "count": int(len(target)),
        "spearman": finite_spearman(target, prediction),
        "pearson": finite_pearson(target, prediction),
    }


def verify_hash_manifest(directory: Path, manifest_name: str):
    manifest = json.loads((directory / manifest_name).read_text())
    checks = {}
    for name, expected in manifest.items():
        path = directory / name
        if not path.exists():
            raise RuntimeError(f"Missing hashed file: {path}")
        observed = sha256_file(path)
        checks[name] = {
            "expected": expected,
            "observed": observed,
            "match": bool(observed == expected),
        }
        if observed != expected:
            raise RuntimeError(f"Hash mismatch: {path}")
    return checks


def evaluate_direction(source, target, pred_dir, target_label_dir):
    pred_file = pred_dir / f"{source.upper()}_TO_{target.upper()}_PREDICTIONS.csv"
    pred = pd.read_csv(pred_file)
    labels = pd.read_csv(target_label_dir / "SOURCE_LABELS.csv")

    if set(labels["dataset"].astype(str).unique()) != {target}:
        raise RuntimeError(
            f"Target label artifact mismatch for {source}->{target}"
        )

    keys = ["dataset", "order", "candidate_id"]
    joined = pred.merge(
        labels,
        on=keys,
        how="inner",
        validate="one_to_one",
    )
    if len(joined) != 384:
        raise RuntimeError(
            f"{source}->{target}: joined count {len(joined)} != 384"
        )

    residual = {}
    reconstructed = {}

    def one_block(frame):
        y_ho = frame["residual_ho"].to_numpy(dtype=float)
        y_true = frame["transformed_target"].to_numpy(dtype=float)
        o2 = frame["exact_o2"].to_numpy(dtype=float)

        r = {}
        q = {
            "EXACT_O2": rank_metrics(y_true, o2),
        }
        for arm in ARMS:
            col = f"pred_{arm.lower()}"
            values = frame[col].to_numpy(dtype=float)
            r[arm] = regression_metrics(y_ho, values)
            q[arm] = rank_metrics(y_true, o2 + values)
        return r, q

    whole_r, whole_q = one_block(joined)
    residual["whole"] = whole_r
    reconstructed["whole"] = whole_q

    for order in (3, 4):
        cell = joined[joined["order"].eq(order)].copy()
        if len(cell) != 192:
            raise RuntimeError(
                f"{source}->{target} order {order}: {len(cell)} != 192"
            )
        rr, qq = one_block(cell)
        residual[f"order_{order}"] = rr
        reconstructed[f"order_{order}"] = qq

    return joined, {
        "source": source,
        "target": target,
        "residual_prediction": residual,
        "reconstructed_fitness": reconstructed,
    }


def gate_direction(block):
    r3 = block["residual_prediction"]["order_3"]
    r4 = block["residual_prediction"]["order_4"]
    q4 = block["reconstructed_fitness"]["order_4"]

    full4 = r4["FULL"]["spearman"]
    low4 = r4["LOW_ONLY"]["spearman"]
    shuffled4 = r4["SHUFFLED_FULL"]["spearman"]
    rep4 = r4["REPRESENTATION_ONLY"]["spearman"]
    full3 = r3["FULL"]["spearman"]
    low3 = r3["LOW_ONLY"]["spearman"]
    full_fit4 = q4["FULL"]["spearman"]
    o2_fit4 = q4["EXACT_O2"]["spearman"]

    values = [
        full4, low4, shuffled4, rep4,
        full3, low3, full_fit4, o2_fit4,
    ]
    if any(v is None for v in values):
        raise RuntimeError("Undefined gate metric.")

    checks = {
        "quad_full_residual_spearman_at_least_0_25": bool(
            full4 >= 0.25
        ),
        "quad_full_beats_low_only_by_0_10": bool(
            full4 >= low4 + 0.10
        ),
        "quad_full_beats_shuffled_by_0_15": bool(
            full4 >= shuffled4 + 0.15
        ),
        "quad_representation_only_at_least_0_15": bool(
            rep4 >= 0.15
        ),
        "quad_full_reconstructed_beats_o2_by_0_10": bool(
            full_fit4 >= o2_fit4 + 0.10
        ),
        "order3_full_not_below_low_only_minus_0_02": bool(
            full3 >= low3 - 0.02
        ),
    }
    checks["all_metric_requirements"] = bool(all(checks.values()))
    return checks


def rejection_condition(block):
    r4 = block["residual_prediction"]["order_4"]
    full4 = r4["FULL"]["spearman"]
    low4 = r4["LOW_ONLY"]["spearman"]
    if full4 is None or low4 is None:
        return False
    return bool(
        full4 < 0.10
        and (full4 - low4) < 0.05
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--representations", required=True)
    parser.add_argument("--gb1-to-trpb", required=True)
    parser.add_argument("--trpb-to-gb1", required=True)
    parser.add_argument("--gb1-labels", required=True)
    parser.add_argument("--trpb-labels", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    rep_dir = Path(args.representations)
    gb_to_tr = Path(args.gb1_to_trpb)
    tr_to_gb = Path(args.trpb_to_gb1)
    gb_labels = Path(args.gb1_labels)
    tr_labels = Path(args.trpb_labels)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    rep_hash_checks = verify_hash_manifest(
        rep_dir, "REPRESENTATION_HASHES.json"
    )
    gb_pred_hash_checks = verify_hash_manifest(
        gb_to_tr, "PREDICTION_HASHES.json"
    )
    tr_pred_hash_checks = verify_hash_manifest(
        tr_to_gb, "PREDICTION_HASHES.json"
    )
    gb_label_hash_checks = verify_hash_manifest(
        gb_labels, "LABEL_HASHES.json"
    )
    tr_label_hash_checks = verify_hash_manifest(
        tr_labels, "LABEL_HASHES.json"
    )

    model_manifest = json.loads(
        (rep_dir / "MODEL_MANIFEST.json").read_text()
    )
    rep_diag = json.loads(
        (rep_dir / "REPRESENTATION_DIAGNOSTICS.json").read_text()
    )
    gb_manifest = json.loads(
        (gb_to_tr / "DIRECTION_MANIFEST.json").read_text()
    )
    tr_manifest = json.loads(
        (tr_to_gb / "DIRECTION_MANIFEST.json").read_text()
    )
    gb_params = json.loads(
        (gb_to_tr / "DECODER_PARAMETERS.json").read_text()
    )
    tr_params = json.loads(
        (tr_to_gb / "DECODER_PARAMETERS.json").read_text()
    )

    weights = model_manifest.get("files", {}).get(
        "model.safetensors", {}
    )

    integrity = {
        "model_revision_exact": bool(
            model_manifest.get("resolved_revision")
            == EXPECTED_MODEL_REVISION
        ),
        "model_weights_hash_exact": bool(
            weights.get("sha256") == EXPECTED_MODEL_SHA256
        ),
        "representation_replay_exact": bool(
            rep_diag.get("deterministic_replay_exact") is True
        ),
        "gb1_to_trpb_direction_exact": bool(
            gb_manifest.get("source_dataset") == "GB1"
            and gb_manifest.get("target_dataset") == "TrpB"
            and gb_manifest.get("target_label_file_present") is False
        ),
        "trpb_to_gb1_direction_exact": bool(
            tr_manifest.get("source_dataset") == "TrpB"
            and tr_manifest.get("target_dataset") == "GB1"
            and tr_manifest.get("target_label_file_present") is False
        ),
        "all_decoder_replays_exact": bool(
            all(
                arm["deterministic_replay_max_abs_difference"] == 0.0
                for params in (gb_params, tr_params)
                for arm in params["arms"].values()
            )
        ),
        "representation_hashes_exact": bool(
            all(x["match"] for x in rep_hash_checks.values())
        ),
        "gb1_to_trpb_hashes_exact": bool(
            all(x["match"] for x in gb_pred_hash_checks.values())
        ),
        "trpb_to_gb1_hashes_exact": bool(
            all(x["match"] for x in tr_pred_hash_checks.values())
        ),
        "gb1_label_hashes_exact": bool(
            all(x["match"] for x in gb_label_hash_checks.values())
        ),
        "trpb_label_hashes_exact": bool(
            all(x["match"] for x in tr_label_hash_checks.values())
        ),
        "nucb_not_consumed": bool(
            rep_diag.get("nucb_consumed") is False
            and gb_manifest.get("nucb_consumed") is False
            and tr_manifest.get("nucb_consumed") is False
        ),
        "phase3_closed": bool(
            rep_diag.get("phase3_opened") is False
            and gb_manifest.get("phase3_opened") is False
            and tr_manifest.get("phase3_opened") is False
        ),
    }
    integrity["all_required_pass"] = bool(all(integrity.values()))

    joined_gb_to_tr, block_gb_to_tr = evaluate_direction(
        "GB1", "TrpB", gb_to_tr, tr_labels
    )
    joined_tr_to_gb, block_tr_to_gb = evaluate_direction(
        "TrpB", "GB1", tr_to_gb, gb_labels
    )

    gate_gb_to_tr = gate_direction(block_gb_to_tr)
    gate_tr_to_gb = gate_direction(block_tr_to_gb)
    gate_gb_to_tr["integrity_all_required_pass"] = integrity[
        "all_required_pass"
    ]
    gate_tr_to_gb["integrity_all_required_pass"] = integrity[
        "all_required_pass"
    ]

    pass_gb_to_tr = bool(all(gate_gb_to_tr.values()))
    pass_tr_to_gb = bool(all(gate_tr_to_gb.values()))
    supported = bool(pass_gb_to_tr and pass_tr_to_gb)

    reject_gb_to_tr = rejection_condition(block_gb_to_tr)
    reject_tr_to_gb = rejection_condition(block_tr_to_gb)

    if supported:
        decision = "FROZEN_EMBEDDING_TRANSFER_SUPPORTED"
        next_action = "PREREGISTER_ONE_ARCHITECTURE_NO_NUCB_YET"
    elif reject_gb_to_tr and reject_tr_to_gb:
        decision = "FROZEN_EMBEDDING_TRANSFER_REJECTED"
        next_action = (
            "DO_NOT_BUILD_FROM_ESM2_8M_HIDDEN_TRANSFER;"
            "NEXT_PRIOR_CLASS_MUST_BE_GENUINELY_DIFFERENT_OR_ADD_MEASUREMENTS"
        )
    else:
        decision = "FROZEN_EMBEDDING_TRANSFER_UNRESOLVED"
        next_action = "DO_NOT_BUILD_ARCHITECTURE"

    metrics = {
        "version": "NABU_FROZEN_EMBEDDING_TRANSFER_METRICS_V1",
        "GB1_to_TrpB": block_gb_to_tr,
        "TrpB_to_GB1": block_tr_to_gb,
        "gates": {
            "GB1_to_TrpB": gate_gb_to_tr,
            "TrpB_to_GB1": gate_tr_to_gb,
            "GB1_to_TrpB_all_required": pass_gb_to_tr,
            "TrpB_to_GB1_all_required": pass_tr_to_gb,
        },
    }
    (out / "TRANSFER_METRICS.json").write_text(
        json.dumps(metrics, indent=2),
        encoding="utf-8",
    )

    decision_doc = {
        "version": "NABU_FROZEN_EMBEDDING_TRANSFER_DECISION_V1",
        "decision": decision,
        "next_action": next_action,
        "root_cause_status": "ROOT_CAUSE_NOT_YET_IDENTIFIED",
        "GB1_to_TrpB_weak_rejection_condition": reject_gb_to_tr,
        "TrpB_to_GB1_weak_rejection_condition": reject_tr_to_gb,
        "integrity_all_required_pass": integrity["all_required_pass"],
        "nucb_consumed": False,
        "phase3_opened": False,
    }
    (out / "DECISION.json").write_text(
        json.dumps(decision_doc, indent=2),
        encoding="utf-8",
    )

    diagnostics = {
        "integrity": integrity,
        "representation_hash_checks": rep_hash_checks,
        "gb1_to_trpb_prediction_hash_checks": gb_pred_hash_checks,
        "trpb_to_gb1_prediction_hash_checks": tr_pred_hash_checks,
        "gb1_label_hash_checks": gb_label_hash_checks,
        "trpb_label_hash_checks": tr_label_hash_checks,
    }
    (out / "EVALUATION_DIAGNOSTICS.json").write_text(
        json.dumps(diagnostics, indent=2),
        encoding="utf-8",
    )

    joined_gb_to_tr = joined_gb_to_tr.copy()
    joined_gb_to_tr.insert(0, "direction_source", "GB1")
    joined_gb_to_tr.insert(1, "direction_target", "TrpB")
    joined_tr_to_gb = joined_tr_to_gb.copy()
    joined_tr_to_gb.insert(0, "direction_source", "TrpB")
    joined_tr_to_gb.insert(1, "direction_target", "GB1")

    pd.concat(
        [joined_gb_to_tr, joined_tr_to_gb],
        ignore_index=True,
    ).to_csv(out / "JOINED_EVIDENCE.csv", index=False)

    hashes = {
        path.name: sha256_file(path)
        for path in sorted(out.iterdir())
        if path.is_file() and path.name != "OUTPUT_HASHES.json"
    }
    (out / "OUTPUT_HASHES.json").write_text(
        json.dumps(hashes, indent=2),
        encoding="utf-8",
    )

    print(json.dumps(decision_doc, indent=2))


if __name__ == "__main__":
    main()
