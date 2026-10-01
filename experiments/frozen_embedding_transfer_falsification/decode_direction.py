from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

SEED = 161
RIDGE_ALPHA = 1.0
LOW_FEATURES = [f"low_feature_{i:02d}" for i in range(15)]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_hash_manifest(directory: Path, manifest_name: str):
    manifest = json.loads((directory / manifest_name).read_text())
    checks = {}
    for name, expected in manifest.items():
        observed = sha256_file(directory / name)
        checks[name] = {
            "expected": expected,
            "observed": observed,
            "match": bool(observed == expected),
        }
        if observed != expected:
            raise RuntimeError(f"Hash mismatch: {directory / name}")
    return checks


def load_feature_table(visible_dir: Path, rep_dir: Path):
    selected = pd.read_csv(visible_dir / "SELECTED_CANDIDATES.csv")
    meta = pd.read_csv(rep_dir / "CANDIDATE_R_HO_METADATA.csv")
    matrix = np.load(
        rep_dir / "CANDIDATE_R_HO.npy",
        allow_pickle=False,
    )

    if len(meta) != len(matrix):
        raise RuntimeError("Representation metadata/matrix mismatch.")

    expected_indices = np.arange(len(meta), dtype=int)
    if not np.array_equal(
        meta["representation_row_index"].to_numpy(dtype=int),
        expected_indices,
    ):
        raise RuntimeError("Representation row index is not canonical.")

    keys = [
        "dataset",
        "order",
        "candidate_id",
        "selection_rank",
        "selection_hash",
    ]
    merged = selected.merge(
        meta[keys + ["representation_row_index"]],
        on=keys,
        how="inner",
        validate="one_to_one",
    )
    if len(merged) != len(selected):
        raise RuntimeError("Selected/representation join mismatch.")

    row_index = merged["representation_row_index"].to_numpy(dtype=int)
    rep = matrix[row_index]

    low = merged[LOW_FEATURES].to_numpy(dtype=float)
    if not np.isfinite(low).all() or not np.isfinite(rep).all():
        raise RuntimeError("Non-finite feature value.")

    return merged, low, rep


def fit_arm(X_source, y_source, X_target, shuffle_labels=False):
    y_fit = np.asarray(y_source, dtype=float).copy()
    if shuffle_labels:
        y_fit = np.random.default_rng(SEED).permutation(y_fit)

    scaler = StandardScaler()
    Xs = scaler.fit_transform(X_source)
    Xt = scaler.transform(X_target)

    model = Ridge(
        alpha=RIDGE_ALPHA,
        solver="svd",
        fit_intercept=True,
    )
    model.fit(Xs, y_fit)
    pred = np.asarray(model.predict(Xt), dtype=float)

    # Exact replay.
    scaler2 = StandardScaler()
    Xs2 = scaler2.fit_transform(X_source)
    Xt2 = scaler2.transform(X_target)
    model2 = Ridge(
        alpha=RIDGE_ALPHA,
        solver="svd",
        fit_intercept=True,
    )
    model2.fit(Xs2, y_fit)
    pred2 = np.asarray(model2.predict(Xt2), dtype=float)

    replay_max_abs = float(np.max(np.abs(pred - pred2)))
    if replay_max_abs != 0.0:
        raise RuntimeError(
            f"Decoder determinism failed: {replay_max_abs}"
        )

    params = {
        "feature_count": int(X_source.shape[1]),
        "scaler_mean": [float(x) for x in scaler.mean_],
        "scaler_scale": [float(x) for x in scaler.scale_],
        "ridge_alpha": RIDGE_ALPHA,
        "ridge_solver": "svd",
        "ridge_fit_intercept": True,
        "ridge_intercept": float(model.intercept_),
        "ridge_coef": [float(x) for x in model.coef_],
        "shuffle_labels": bool(shuffle_labels),
        "deterministic_replay_max_abs_difference": replay_max_abs,
    }
    return pred, params


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--visible", required=True)
    parser.add_argument("--representations", required=True)
    parser.add_argument("--source-labels", required=True)
    parser.add_argument("--source-dataset", required=True, choices=["GB1", "TrpB"])
    parser.add_argument("--target-dataset", required=True, choices=["GB1", "TrpB"])
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    if args.source_dataset == args.target_dataset:
        raise RuntimeError("Source and target datasets must differ.")

    visible_dir = Path(args.visible)
    rep_dir = Path(args.representations)
    label_dir = Path(args.source_labels)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    visible_hash_checks = verify_hash_manifest(
        visible_dir,
        "VISIBLE_HASHES.json",
    )
    rep_hash_checks = verify_hash_manifest(
        rep_dir,
        "REPRESENTATION_HASHES.json",
    )
    label_hash_checks = verify_hash_manifest(
        label_dir,
        "LABEL_HASHES.json",
    )

    table, low, rep = load_feature_table(
        visible_dir,
        rep_dir,
    )
    labels = pd.read_csv(label_dir / "SOURCE_LABELS.csv")

    if set(labels["dataset"].astype(str).unique()) != {args.source_dataset}:
        raise RuntimeError(
            "Source-label artifact contains unexpected dataset."
        )
    if args.target_dataset in set(labels["dataset"].astype(str)):
        raise RuntimeError("Target labels visible to prediction job.")

    keys = ["dataset", "order", "candidate_id"]
    source_mask = table["dataset"].astype(str).eq(args.source_dataset)
    target_mask = table["dataset"].astype(str).eq(args.target_dataset)

    source_table = table[source_mask].copy()
    target_table = table[target_mask].copy()
    source_low = low[source_mask.to_numpy()]
    target_low = low[target_mask.to_numpy()]
    source_rep = rep[source_mask.to_numpy()]
    target_rep = rep[target_mask.to_numpy()]

    source_joined = source_table.merge(
        labels,
        on=keys,
        how="inner",
        validate="one_to_one",
    )
    if len(source_joined) != len(source_table):
        raise RuntimeError("Source labels do not cover source rows exactly.")

    # Preserve the source feature order from source_table.
    source_index = {
        (str(r.dataset), int(r.order), str(r.candidate_id)): i
        for i, r in enumerate(source_table.itertuples(index=False))
    }
    source_order = np.asarray(
        [
            source_index[
                (str(r.dataset), int(r.order), str(r.candidate_id))
            ]
            for r in source_joined.itertuples(index=False)
        ],
        dtype=int,
    )
    source_low = source_low[source_order]
    source_rep = source_rep[source_order]
    y_source = source_joined["residual_ho"].to_numpy(dtype=float)

    if len(source_joined) != 384 or len(target_table) != 384:
        raise RuntimeError(
            f"Direction row contract failed: source={len(source_joined)} "
            f"target={len(target_table)}"
        )

    full_source = np.hstack([source_low, source_rep])
    full_target = np.hstack([target_low, target_rep])

    arms = {}
    predictions = {}

    predictions["LOW_ONLY"], arms["LOW_ONLY"] = fit_arm(
        source_low,
        y_source,
        target_low,
        shuffle_labels=False,
    )
    predictions["REPRESENTATION_ONLY"], arms["REPRESENTATION_ONLY"] = fit_arm(
        source_rep,
        y_source,
        target_rep,
        shuffle_labels=False,
    )
    predictions["FULL"], arms["FULL"] = fit_arm(
        full_source,
        y_source,
        full_target,
        shuffle_labels=False,
    )
    predictions["SHUFFLED_FULL"], arms["SHUFFLED_FULL"] = fit_arm(
        full_source,
        y_source,
        full_target,
        shuffle_labels=True,
    )

    pred_df = target_table[
        [
            "dataset",
            "order",
            "candidate_id",
            "selection_rank",
            "selection_hash",
        ]
    ].copy()
    for name, values in predictions.items():
        pred_df[f"pred_{name.lower()}"] = values

    pred_path = out / (
        f"{args.source_dataset.upper()}_TO_"
        f"{args.target_dataset.upper()}_PREDICTIONS.csv"
    )
    pred_df.to_csv(pred_path, index=False)

    parameters = {
        "version": "NABU_FROZEN_EMBEDDING_DECODER_V1",
        "source_dataset": args.source_dataset,
        "target_dataset": args.target_dataset,
        "source_rows": int(len(source_joined)),
        "target_rows": int(len(target_table)),
        "source_order_counts": {
            str(k): int(v)
            for k, v in source_joined.groupby("order").size().items()
        },
        "target_order_counts": {
            str(k): int(v)
            for k, v in target_table.groupby("order").size().items()
        },
        "low_feature_count": int(source_low.shape[1]),
        "representation_feature_count": int(source_rep.shape[1]),
        "full_feature_count": int(full_source.shape[1]),
        "arms": arms,
        "nucb_consumed": False,
        "phase3_opened": False,
    }
    (out / "DECODER_PARAMETERS.json").write_text(
        json.dumps(parameters, indent=2),
        encoding="utf-8",
    )

    manifest = {
        "version": "NABU_FROZEN_EMBEDDING_DIRECTION_FREEZE_V1",
        "source_dataset": args.source_dataset,
        "target_dataset": args.target_dataset,
        "seed": SEED,
        "ridge_alpha": RIDGE_ALPHA,
        "ridge_solver": "svd",
        "source_label_file_sha256": sha256_file(
            label_dir / "SOURCE_LABELS.csv"
        ),
        "target_label_file_present": False,
        "visible_hash_checks": visible_hash_checks,
        "representation_hash_checks": rep_hash_checks,
        "source_label_hash_checks": label_hash_checks,
        "prediction_file": pred_path.name,
        "nucb_consumed": False,
        "phase3_opened": False,
    }
    (out / "DIRECTION_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )

    hashes = {
        path.name: sha256_file(path)
        for path in sorted(out.iterdir())
        if path.is_file() and path.name != "PREDICTION_HASHES.json"
    }
    (out / "PREDICTION_HASHES.json").write_text(
        json.dumps(hashes, indent=2),
        encoding="utf-8",
    )

    print(json.dumps({
        "source": args.source_dataset,
        "target": args.target_dataset,
        "source_rows": int(len(source_joined)),
        "target_rows": int(len(target_table)),
        "prediction_file": pred_path.name,
        "all_decoder_replays_exact": bool(
            all(
                arm["deterministic_replay_max_abs_difference"] == 0.0
                for arm in arms.values()
            )
        ),
        "nucb_consumed": False,
    }, indent=2))


if __name__ == "__main__":
    main()
