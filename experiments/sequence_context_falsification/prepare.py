from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

SEED = 161
EPS = 1e-12
SAMPLE_N = 192
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

TRPB_EXPECTED_MD5 = "a611408d2db06907a023ddc8a1d94c12"
TRPB_EXPECTED_COUNTS = {"fit": 8633, "validation": 2158, "test": 217507}

GB1_WT = "VDGV"
GB1_POSITIONS = (39, 40, 41, 54)
GB1_REFERENCE = "MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE"
GB1_EXPECTED_SIZE = 2993024
GB1_EXPECTED_GIT_BLOB = "a1e9f5146a9441210eaf88accfa660f02d80cbb7"


def hash_file(path: Path, algo: str) -> str:
    h = hashlib.new(algo)
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def git_blob_sha(path: Path) -> str:
    data = path.read_bytes()
    payload = f"blob {len(data)}\0".encode("utf-8") + data
    return hashlib.sha1(payload).hexdigest()


def validation_mask(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series.dtype):
        return series.fillna(False).astype(bool)
    text = series.fillna("").astype(str).str.strip().str.lower()
    return text.isin({"true", "1", "yes", "y", "t"})


def derive_reference(sequences: list[str]) -> str:
    lengths = {len(x) for x in sequences}
    if len(lengths) != 1:
        raise RuntimeError("TrpB fit sequences are not fixed length.")
    if any(not set(x).issubset(STANDARD_AA) for x in sequences):
        raise RuntimeError("TrpB contains non-standard amino acids.")
    length = next(iter(lengths))
    return "".join(
        Counter(x[i] for x in sequences).most_common(1)[0][0]
        for i in range(length)
    )


def mutation_set(sequence: str, reference: str) -> tuple[str, ...]:
    sequence = str(sequence).strip().upper()
    if len(sequence) != len(reference):
        raise RuntimeError("Sequence/reference length mismatch.")
    return tuple(
        f"{source}{i + 1}{target}"
        for i, (source, target) in enumerate(zip(reference, sequence))
        if source != target
    )


def apply_mutations(reference: str, muts: tuple[str, ...]) -> str:
    chars = list(reference)
    seen = set()
    for token in muts:
        source = token[0]
        target = token[-1]
        pos = int(token[1:-1]) - 1
        if pos in seen:
            raise RuntimeError(f"Duplicate mutation position: {muts}")
        seen.add(pos)
        if chars[pos] != source:
            raise RuntimeError(
                f"Mutation source mismatch at {pos + 1}: "
                f"{chars[pos]} != {source}"
            )
        chars[pos] = target
    return "".join(chars)


class RobustAsinhTransform:
    def fit(self, y):
        y = np.asarray(y, dtype=float)
        self.center = float(np.median(y))
        mad = float(np.median(np.abs(y - self.center)))
        self.scale_source = "mad"
        if mad <= EPS:
            mad = float(np.std(y))
            self.scale_source = "std_fallback"
        if mad <= EPS:
            mad = 1.0
            self.scale_source = "unit_fallback"
        self.scale = mad
        return self

    def transform(self, y):
        y = np.asarray(y, dtype=float)
        return np.arcsinh((y - self.center) / self.scale)


def unique_lookup(mutation_sets, values):
    out = {}
    for muts, value in zip(mutation_sets, values):
        key = tuple(muts)
        if key in out:
            raise RuntimeError(f"Duplicate fit mutation set: {key}")
        out[key] = float(value)
    return out


def exact_components(lookup):
    if () not in lookup:
        raise RuntimeError("WT missing from low-order fit.")
    wt = float(lookup[()])
    e1 = {
        muts[0]: float(value - wt)
        for muts, value in lookup.items()
        if len(muts) == 1
    }
    e2 = {}
    for muts, value in lookup.items():
        if len(muts) != 2:
            continue
        left, right = muts
        if left in e1 and right in e1:
            e2[muts] = float(value - wt - e1[left] - e1[right])
    return wt, e1, e2


def low_order_vector(muts, wt, e1, e2, sigma_low):
    if not all(m in e1 for m in muts):
        return None
    pairs = [tuple(x) for x in combinations(muts, 2)]
    if not all(p in e2 for p in pairs):
        return None

    e1v = np.asarray([e1[m] for m in muts], dtype=float)
    e2v = np.asarray([e2[p] for p in pairs], dtype=float)
    o1 = float(wt + np.sum(e1v))
    o2 = float(o1 + np.sum(e2v))
    features = [
        o1 / sigma_low,
        o2 / sigma_low,
        float(np.sum(e1v)) / sigma_low,
        float(np.mean(np.abs(e1v))) / sigma_low,
        float(np.max(np.abs(e1v))) / sigma_low,
        float(np.std(e1v)) / sigma_low,
        float(np.min(e1v)) / sigma_low,
        float(np.max(e1v)) / sigma_low,
        float(np.sum(e2v)) / sigma_low,
        float(np.mean(np.abs(e2v))) / sigma_low,
        float(np.max(np.abs(e2v))) / sigma_low,
        float(np.std(e2v)) / sigma_low,
        float(np.min(e2v)) / sigma_low,
        float(np.max(e2v)) / sigma_low,
        float(np.mean(e2v > 0.0)),
    ]
    return o2, features


def selection_hash(dataset: str, order: int, candidate_id: str) -> str:
    payload = (
        f"NABU_SEQ_CONTEXT_V1|{SEED}|{dataset}|{order}|{candidate_id}"
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def load_trpb(path: Path):
    if hash_file(path, "md5") != TRPB_EXPECTED_MD5:
        raise RuntimeError("TrpB source MD5 mismatch.")

    frame = pd.read_csv(path, compression="gzip")
    required = {"sequence", "target", "set", "validation"}
    if required - set(frame.columns):
        raise RuntimeError("TrpB source columns missing.")

    frame = frame.copy()
    frame["sequence"] = frame["sequence"].astype(str).str.strip().str.upper()
    frame["target"] = pd.to_numeric(frame["target"], errors="raise")

    set_text = frame["set"].astype(str).str.strip().str.lower()
    val = validation_mask(frame["validation"])
    fit = frame[set_text.eq("train") & ~val].copy()
    validation = frame[set_text.eq("train") & val].copy()
    test = frame[set_text.eq("test")].copy()

    counts = {
        "fit": int(len(fit)),
        "validation": int(len(validation)),
        "test": int(len(test)),
    }
    if counts != TRPB_EXPECTED_COUNTS:
        raise RuntimeError(f"TrpB split drift: {counts}")

    reference = derive_reference(fit["sequence"].tolist())
    fit_sets = [mutation_set(x, reference) for x in fit["sequence"]]
    test_sets = [mutation_set(x, reference) for x in test["sequence"]]

    return {
        "name": "TrpB",
        "fit_ids": fit["sequence"].tolist(),
        "fit_sequences": fit["sequence"].tolist(),
        "fit_sets": fit_sets,
        "fit_y": fit["target"].to_numpy(dtype=float),
        "test_ids": test["sequence"].tolist(),
        "test_sequences": test["sequence"].tolist(),
        "test_sets": test_sets,
        "test_y": test["target"].to_numpy(dtype=float),
        "reference": reference,
        "source": {
            "md5": hash_file(path, "md5"),
            "sha256": hash_file(path, "sha256"),
            "counts": counts,
        },
    }


def gb1_mutations(genotype: str) -> tuple[str, ...]:
    genotype = str(genotype).strip().upper()
    if len(genotype) != 4:
        raise RuntimeError(f"Invalid GB1 genotype: {genotype}")
    return tuple(
        f"{source}{pos}{target}"
        for source, target, pos in zip(
            GB1_WT, genotype, GB1_POSITIONS
        )
        if source != target
    )


def gb1_full_sequence(genotype: str) -> str:
    muts = gb1_mutations(genotype)
    return apply_mutations(GB1_REFERENCE, muts)


def load_gb1(path: Path):
    if path.stat().st_size != GB1_EXPECTED_SIZE:
        raise RuntimeError("GB1 source size drift.")
    if git_blob_sha(path) != GB1_EXPECTED_GIT_BLOB:
        raise RuntimeError("GB1 source git-blob identity drift.")

    wt_from_sequence = "".join(
        GB1_REFERENCE[p - 1] for p in GB1_POSITIONS
    )
    if wt_from_sequence != GB1_WT:
        raise RuntimeError(
            f"GB1 reference invariant failed: {wt_from_sequence}"
        )

    frame = pd.read_csv(path)
    frame = frame[["variant", "fitness"]].dropna().copy()
    frame["variant"] = frame["variant"].astype(str).str.strip().str.upper()
    frame["fitness"] = pd.to_numeric(frame["fitness"], errors="raise")
    if frame["variant"].duplicated().any():
        raise RuntimeError("GB1 duplicate variants.")

    frame["mutation_set"] = frame["variant"].map(gb1_mutations)
    frame["order"] = frame["mutation_set"].map(len)
    fit = frame[frame["order"].le(2)].copy()
    test = frame[frame["order"].ge(3)].copy()

    return {
        "name": "GB1",
        "fit_ids": fit["variant"].tolist(),
        "fit_sequences": fit["variant"].map(gb1_full_sequence).tolist(),
        "fit_sets": fit["mutation_set"].tolist(),
        "fit_y": fit["fitness"].to_numpy(dtype=float),
        "test_ids": test["variant"].tolist(),
        "test_sequences": test["variant"].map(gb1_full_sequence).tolist(),
        "test_sets": test["mutation_set"].tolist(),
        "test_y": test["fitness"].to_numpy(dtype=float),
        "reference": GB1_REFERENCE,
        "source": {
            "size": int(path.stat().st_size),
            "git_blob_sha": git_blob_sha(path),
            "sha256": hash_file(path, "sha256"),
            "row_count": int(len(frame)),
        },
    }


def prepare_dataset(dataset, visible_rows, reveal_rows, diagnostics):
    transform = RobustAsinhTransform().fit(dataset["fit_y"])
    fit_y = transform.transform(dataset["fit_y"])
    test_y = transform.transform(dataset["test_y"])

    lookup = unique_lookup(dataset["fit_sets"], fit_y)
    wt, e1, e2 = exact_components(lookup)
    sigma_low = float(np.std(fit_y))
    if sigma_low <= EPS:
        raise RuntimeError(f"{dataset['name']}: zero low-order scale.")

    eligible = {3: [], 4: []}
    for idx, (cid, seq, muts) in enumerate(
        zip(
            dataset["test_ids"],
            dataset["test_sequences"],
            dataset["test_sets"],
        )
    ):
        order = len(muts)
        if order not in eligible:
            continue
        built = low_order_vector(muts, wt, e1, e2, sigma_low)
        if built is None:
            continue
        o2, features = built
        eligible[order].append(
            {
                "idx": idx,
                "candidate_id": str(cid),
                "sequence": seq,
                "muts": tuple(muts),
                "o2": float(o2),
                "features": features,
                "selection_hash": selection_hash(
                    dataset["name"], order, str(cid)
                ),
            }
        )

    for order in (3, 4):
        rows = sorted(
            eligible[order], key=lambda x: x["selection_hash"]
        )
        if len(rows) < SAMPLE_N:
            raise RuntimeError(
                f"{dataset['name']} order {order}: "
                f"only {len(rows)} eligible < {SAMPLE_N}"
            )
        selected = rows[:SAMPLE_N]
        for rank, row in enumerate(selected):
            idx = row["idx"]
            base = {
                "dataset": dataset["name"],
                "order": order,
                "selection_rank": rank,
                "candidate_id": row["candidate_id"],
                "selection_hash": row["selection_hash"],
                "reference_sequence": dataset["reference"],
                "candidate_sequence": row["sequence"],
                "mutation_tokens": json.dumps(list(row["muts"])),
                "sigma_low": sigma_low,
            }
            for j, value in enumerate(row["features"]):
                base[f"low_feature_{j:02d}"] = float(value)
            visible_rows.append(base)

            reveal_rows.append(
                {
                    "dataset": dataset["name"],
                    "order": order,
                    "candidate_id": row["candidate_id"],
                    "raw_target": float(dataset["test_y"][idx]),
                    "transformed_target": float(test_y[idx]),
                    "exact_o2": row["o2"],
                    "residual_ho": float(test_y[idx] - row["o2"]),
                }
            )

    diagnostics[dataset["name"]] = {
        "reference_length": int(len(dataset["reference"])),
        "reference_sha256": hashlib.sha256(
            dataset["reference"].encode("utf-8")
        ).hexdigest(),
        "fit_count": int(len(dataset["fit_y"])),
        "test_count": int(len(dataset["test_y"])),
        "fit_order_counts": dict(
            Counter(str(len(x)) for x in dataset["fit_sets"])
        ),
        "test_order_counts": dict(
            Counter(str(len(x)) for x in dataset["test_sets"])
        ),
        "eligible_counts": {
            str(order): int(len(eligible[order]))
            for order in (3, 4)
        },
        "selected_counts": {
            "3": SAMPLE_N,
            "4": SAMPLE_N,
        },
        "low_order_components": {
            "e1_count": int(len(e1)),
            "e2_count": int(len(e2)),
            "sigma_low": sigma_low,
            "robust_asinh_center": float(transform.center),
            "robust_asinh_scale": float(transform.scale),
            "robust_asinh_scale_source": transform.scale_source,
        },
        "source": dataset["source"],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--trpb", required=True)
    parser.add_argument("--gb1", required=True)
    parser.add_argument("--visible-out", required=True)
    parser.add_argument("--reveal-out", required=True)
    args = parser.parse_args()

    visible_out = Path(args.visible_out)
    reveal_out = Path(args.reveal_out)
    visible_out.mkdir(parents=True, exist_ok=True)
    reveal_out.mkdir(parents=True, exist_ok=True)

    trpb = load_trpb(Path(args.trpb))
    gb1 = load_gb1(Path(args.gb1))

    visible_rows = []
    reveal_rows = []
    diagnostics = {}

    prepare_dataset(trpb, visible_rows, reveal_rows, diagnostics)
    prepare_dataset(gb1, visible_rows, reveal_rows, diagnostics)

    visible = pd.DataFrame(visible_rows)
    reveal = pd.DataFrame(reveal_rows)

    expected_total = 2 * 2 * SAMPLE_N
    if len(visible) != expected_total or len(reveal) != expected_total:
        raise RuntimeError("Selected/reveal row count contract failed.")

    visible.to_csv(
        visible_out / "SELECTED_CANDIDATES.csv",
        index=False,
    )
    reveal.to_csv(
        reveal_out / "REVEAL_TARGETS.csv",
        index=False,
    )

    (visible_out / "LOW_ORDER_DIAGNOSTICS.json").write_text(
        json.dumps(diagnostics, indent=2),
        encoding="utf-8",
    )

    source_manifest = {
        "version": "NABU_SEQUENCE_CONTEXT_FALSIFICATION_PREP_V1",
        "seed": SEED,
        "sample_n_per_cell": SAMPLE_N,
        "selection_rule": (
            "sha256(NABU_SEQ_CONTEXT_V1|161|dataset|order|candidate_id)"
        ),
        "trpb": trpb["source"],
        "gb1": gb1["source"],
        "gb1_reference": {
            "pdb": "3GB1",
            "sequence": GB1_REFERENCE,
            "positions": list(GB1_POSITIONS),
            "wt_genotype": GB1_WT,
        },
        "nucb_consumed": False,
        "phase3_opened": False,
    }
    (visible_out / "SOURCE_MANIFEST.json").write_text(
        json.dumps(source_manifest, indent=2),
        encoding="utf-8",
    )

    visible_hashes = {
        p.name: hash_file(p, "sha256")
        for p in sorted(visible_out.iterdir())
        if p.is_file()
    }
    reveal_hashes = {
        p.name: hash_file(p, "sha256")
        for p in sorted(reveal_out.iterdir())
        if p.is_file()
    }
    (visible_out / "VISIBLE_HASHES.json").write_text(
        json.dumps(visible_hashes, indent=2),
        encoding="utf-8",
    )
    (reveal_out / "REVEAL_HASHES.json").write_text(
        json.dumps(reveal_hashes, indent=2),
        encoding="utf-8",
    )

    print(json.dumps({
        "selected_rows": int(len(visible)),
        "cells": visible.groupby(["dataset", "order"]).size().to_dict(),
        "nucb_consumed": False,
    }, indent=2, default=str))


if __name__ == "__main__":
    main()
