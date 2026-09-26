from __future__ import annotations

import hashlib

import numpy as np
from scipy.stats import rankdata


def rank_hash(candidate_ids, prediction) -> str:
    candidate_ids = np.asarray(candidate_ids, dtype=str)
    prediction = np.asarray(prediction, dtype=float)
    order = np.lexsort((candidate_ids, -prediction))
    payload = "\n".join(candidate_ids[order].tolist())
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


class _FenwickMax:
    def __init__(self, size: int):
        self.tree = np.zeros(size + 1, dtype=np.int64)

    def update(self, index: int, value: int) -> None:
        index += 1
        size = len(self.tree)
        while index < size:
            if value > self.tree[index]:
                self.tree[index] = value
            index += index & -index

    def query_prefix(self, index: int) -> int:
        index += 1
        result = 0
        while index > 0:
            if self.tree[index] > result:
                result = int(self.tree[index])
            index -= index & -index
        return result


def pareto_consensus_score_scalable(
    b2_score: np.ndarray,
    context_score: np.ndarray,
    candidate_ids,
):
    b2_score = np.asarray(b2_score, dtype=float)
    context_score = np.asarray(context_score, dtype=float)
    candidate_ids = np.asarray(candidate_ids, dtype=str)

    n = len(b2_score)
    if n == 0:
        raise ValueError("Consensus requires at least one candidate.")
    if len(context_score) != n or len(candidate_ids) != n:
        raise ValueError("Consensus input length mismatch.")
    if not np.isfinite(b2_score).all() or not np.isfinite(context_score).all():
        raise ValueError("Consensus scores must be finite.")

    b2_rank = rankdata(-b2_score, method="average") / float(n)
    context_rank = rankdata(-context_score, method="average") / float(n)

    context_values = np.unique(context_rank)
    context_position = {
        float(value): index
        for index, value in enumerate(context_values.tolist())
    }

    order = np.lexsort((context_rank, b2_rank))
    layer = np.full(n, -1, dtype=np.int64)
    fenwick = _FenwickMax(len(context_values))

    start = 0
    while start < n:
        first = order[start]
        current_b2 = b2_rank[first]
        current_context = context_rank[first]

        end = start + 1
        while end < n:
            row = order[end]
            if (
                b2_rank[row] != current_b2
                or context_rank[row] != current_context
            ):
                break
            end += 1

        position = context_position[float(current_context)]
        current_layer = fenwick.query_prefix(position)
        group = order[start:end]
        layer[group] = current_layer

        fenwick.update(position, current_layer + 1)
        start = end

    if np.any(layer < 0):
        raise RuntimeError("Failed to assign every Pareto layer.")

    mean_rank = 0.5 * (b2_rank + context_rank)
    final_order = np.lexsort((candidate_ids, mean_rank, layer))
    score = np.empty(n, dtype=float)
    score[final_order] = np.arange(n, 0, -1, dtype=float)

    return score, layer, b2_rank, context_rank
