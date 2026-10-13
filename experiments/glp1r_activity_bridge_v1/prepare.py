"""Audit source records and separate external identities from external labels."""
from __future__ import annotations

from collections import Counter

import numpy as np
import pandas as pd

from common import HERE, ALPHABET, folds_for, read_json, sequence_groups, sha256, write_json


def decompose(raw, aligned):
    core = aligned.replace("-", "")
    if len(aligned) != 30 or raw.count(core) != 1:
        raise ValueError("Ambiguous or noncontiguous source alignment")
    start = raw.index(core)
    prefix, suffix = raw[:start], raw[start + len(core):]
    assert prefix + core + suffix == raw
    return prefix, suffix


def prepare():
    config = read_json(HERE / "protocol.json")
    raw_dir, out = HERE / "data/raw", HERE / "data/prepared"
    if (HERE / "results/PREDICTION_FREEZE.json").exists():
        raise RuntimeError("Preparation cannot overwrite an outcome-frozen run")
    for name, expected in config["source_sha256"].items():
        assert sha256(raw_dir / name) == expected, name
    train = pd.read_excel(raw_dir / "training_data.xlsx", sheet_name="dataset")
    aligned = pd.read_excel(raw_dir / "training_data.xlsx", sheet_name="alignment")
    external = pd.read_excel(raw_dir / "reference_data.xlsx", sheet_name="2009Day")
    assert len(train) == len(aligned) == 125 and len(external) == 19
    assert train.pep_ID.tolist() == aligned.pep_ID.tolist()
    assert not train.sequence.duplicated().any()
    assert train.sequence.str.len().eq(train.length).all()
    parts = [decompose(s, a) for s, a in zip(train.sequence, aligned.sequence)]
    prefix_width = max(len(p) for p, _ in parts)
    suffix_width = max(len(s) for _, s in parts)

    def encode(raw, core):
        prefix, suffix = decompose(raw, core)
        if len(prefix) > prefix_width or len(suffix) > suffix_width:
            raise ValueError("External sequence exceeds frozen feature schema")
        encoded = prefix.rjust(prefix_width, "-") + core + suffix.ljust(suffix_width, "-")
        assert encoded.replace("-", "") == raw
        assert set(encoded) <= set(ALPHABET)
        return encoded, len(prefix), len(core.replace("-", "")), len(suffix)

    rows = []
    for t, a in zip(train.to_dict("records"), aligned.to_dict("records")):
        enc, n, c, s = encode(t["sequence"], a["sequence"])
        rows.append({"candidate_id": t["pep_ID"], "sequence": t["sequence"],
                     "encoded_sequence": enc, "prefix_length": n, "core_length": c,
                     "suffix_length": s, "target_pEC50": -float(t["EC50_LOG_T2"])})
    dev = pd.DataFrame(rows)
    identities, truth = [], []
    for record in external.to_dict("records"):
        raw = record["sequence"].replace("-", "")
        enc, n, c, s = encode(raw, record["sequence"])
        identities.append({"candidate_id": record["pep_ID"], "sequence": raw,
                           "encoded_sequence": enc, "prefix_length": n,
                           "core_length": c, "suffix_length": s})
        truth.append({"candidate_id": record["pep_ID"], "sequence": raw,
                      "target_pEC50": -float(record["ec_50_log_GLP1R"])})
    ext = pd.DataFrame(identities)
    assert not set(dev.sequence).intersection(ext.sequence)
    groups, distances = sequence_groups(dev.sequence.tolist(), config["cluster_edit_distance"])
    dev["sequence_group"] = groups
    dev["outer_fold"] = -1
    minimum_cross_distances = []
    for fold, (fit, hold) in enumerate(folds_for(groups, config["outer_folds"])):
        dev.loc[hold, "outer_fold"] = fold
        minimum = int(distances[np.ix_(fit, hold)].min())
        assert minimum > config["cluster_edit_distance"]
        minimum_cross_distances.append(minimum)
    assert (dev.outer_fold >= 0).all()
    unit_error = {}
    for suffix in ["T1", "T2"]:
        error = np.abs(np.log10(train["EC50_" + suffix] * 1e-12) - train["EC50_LOG_" + suffix])
        unit_error[suffix] = float(error.max())
        assert error.max() < 0.01
    out.mkdir(parents=True, exist_ok=True)
    dev.to_csv(out / "development.csv", index=False)
    ext.to_csv(out / "external_identities.csv", index=False)
    pd.DataFrame(truth).to_csv(out / "external_truth.csv", index=False)
    dev[["candidate_id", "sequence_group", "outer_fold"]].to_csv(out / "split_manifest.csv", index=False)
    audit = {
        "source_commit": config["source_commit"], "source_sha256": config["source_sha256"],
        "development_records": len(dev), "external_records": len(ext),
        "external_unique_sequences": int(ext.sequence.nunique()),
        "source_alignment_width": 30, "lossless_feature_width": int(dev.encoded_sequence.str.len().iloc[0]),
        "prefix_width": prefix_width, "suffix_width": suffix_width,
        "development_length_counts": {str(k): v for k, v in sorted(Counter(train.sequence.str.len()).items())},
        "source_alignment_omits_residues_records": int(sum(bool(p or s) for p, s in parts)),
        "source_alignment_sequence_collisions": int(aligned.sequence.duplicated().sum()),
        "lossless_sequence_collisions": int(dev.encoded_sequence.duplicated().sum()),
        "omitted_prefix_records": [train.iloc[i].pep_ID for i, (p, _) in enumerate(parts) if p],
        "all_sequences_reconstructed_exactly": True, "all_rows_retained": True,
        "development_external_exact_overlaps": 0,
        "external_sequence_equivalent_record_ids": ext.loc[ext.sequence.duplicated(False), "candidate_id"].tolist(),
        "sequence_component_sizes": sorted(Counter(groups).values(), reverse=True),
        "minimum_outer_fold_train_test_edit_distances": minimum_cross_distances,
        "training_raw_pM_log_molar_max_rounding_error": unit_error,
        "external_labels_not_used_for_model_selection": True,
        "endpoint_limitation": "Reported EC50 only; Emax and complete chemical modification metadata unavailable",
        "prepared_sha256": {p.name: sha256(p) for p in sorted(out.glob("*.csv"))},
    }
    write_json(HERE / "DATA_AUDIT.json", audit)
    print({k: audit[k] for k in ["development_records", "external_records", "external_unique_sequences", "source_alignment_omits_residues_records", "lossless_feature_width"]})


if __name__ == "__main__":
    prepare()
