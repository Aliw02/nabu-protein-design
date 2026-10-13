"""Shared deterministic data and evidence helpers for the GLP1R audit."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from sklearn.model_selection import GroupKFold

HERE = Path(__file__).resolve().parent
ALPHABET = "-ACDEFGHIKLMNPQRSTVWY"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def read_json(path):
    return json.loads(Path(path).read_text())


def levenshtein(left, right):
    previous = list(range(len(right) + 1))
    for i, a in enumerate(left, 1):
        current = [i]
        for j, b in enumerate(right, 1):
            current.append(min(current[-1] + 1, previous[j] + 1, previous[j - 1] + (a != b)))
        previous = current
    return previous[-1]


def sequence_groups(sequences, threshold):
    parents = list(range(len(sequences)))

    def root(index):
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    distances = np.zeros((len(sequences), len(sequences)), dtype=int)
    for i, left in enumerate(sequences):
        for j in range(i):
            distance = levenshtein(left, sequences[j])
            distances[i, j] = distances[j, i] = distance
            if distance <= threshold:
                parents[root(i)] = root(j)
    members = {}
    for i in range(len(sequences)):
        members.setdefault(root(i), []).append(i)
    ordered = sorted(members.values(), key=lambda ids: min(sequences[i] for i in ids))
    groups = np.zeros(len(sequences), dtype=int)
    for group, ids in enumerate(ordered):
        groups[ids] = group
    return groups, distances


def folds_for(groups, n_splits):
    groups = np.asarray(groups)
    if len(np.unique(groups)) < n_splits:
        raise ValueError("Insufficient independent sequence components for frozen fold count")
    return list(GroupKFold(n_splits=n_splits).split(np.zeros(len(groups)), groups=groups))


def codes_from(frame):
    sequences = frame["encoded_sequence"].tolist()
    if len(set(map(len, sequences))) != 1:
        raise ValueError("Feature widths differ")
    return np.array([[ALPHABET.index(c) for c in s] for s in sequences], dtype=np.int16)


def source_manifest():
    paths = sorted(HERE.glob("*.py")) + [HERE / "PROTOCOL.md", HERE / "protocol.json", HERE / "requirements.txt"]
    legacy = HERE.parent / "v8_1_crossfit_higher_order_dev/run_v8_1_crossfit_dev.py"
    return {**{str(p.relative_to(HERE)): sha256(p) for p in paths},
            "../v8_1_crossfit_higher_order_dev/run_v8_1_crossfit_dev.py": sha256(legacy)}
