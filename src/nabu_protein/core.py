"""
NABU Protein Core Architecture
==============================
Mathematical implementation of Structured State Memory, Pairwise Epistatic Residuals,
and Empirical Bayesian Shrinkage for Combinatorial Protein Fitness Landscapes.
"""

from __future__ import annotations

import re
from itertools import combinations
from typing import Dict, Iterable, List, Literal, Optional, Sequence, Set, Tuple, Union

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def fnv1a_hash(s: str) -> int:
    """Computes deterministic 32-bit FNV-1a hash of a string."""
    h = 2166136261
    for b in s.encode("utf-8"):
        h = ((h ^ b) * 16777619) & 0xFFFFFFFF
    return h


def parse_mutations(mutant_str: str) -> Tuple[str, ...]:
    """
    Parses amino acid substitutions from a mutation string (e.g. 'A10C:B20D' or 'A10C/B20D').
    Returns a sorted tuple of mutations ordered by residue position.
    """
    matches = re.findall(r"[A-Z]\d+[A-Z]", str(mutant_str))
    return tuple(sorted(matches, key=lambda x: int(x[1:-1])))


def percentile_rank(population: Sequence[float], value: float) -> float:
    """Calculates empirical percentile rank of a value within a distribution."""
    arr = np.asarray(population)
    return float((np.sum(arr < value) + 0.5 * np.sum(arr == value)) / len(arr))


class ComponentMemoryModel:
    """
    Additive structured state memory baseline.
    Estimates global mean and empirical Bayesian shrunk single-mutation additive effects.
    """

    def __init__(self, shrinkage_prior: float = 1.0):
        self.shrinkage_prior = shrinkage_prior
        self.global_mean: float = 0.0
        self.component_effects: Dict[str, float] = {}

    def fit(self, mutation_sets: Sequence[Tuple[str, ...]], fitness_scores: Sequence[float]) -> ComponentMemoryModel:
        scores = np.asarray(fitness_scores, dtype=float)
        self.global_mean = float(np.mean(scores))
        
        all_mutations: Set[str] = {m for ms in mutation_sets for m in ms}
        self.component_effects = {}
        
        for m in all_mutations:
            indices = [i for i, ms in enumerate(mutation_sets) if m in ms]
            n = len(indices)
            if n > 0:
                subset_mean = float(np.mean(scores[indices]))
                shrinkage = n / (n + self.shrinkage_prior)
                self.component_effects[m] = (subset_mean - self.global_mean) * shrinkage

        return self

    def predict_one(self, mutation_set: Tuple[str, ...]) -> float:
        return self.global_mean + sum(self.component_effects.get(m, 0.0) for m in mutation_set)

    def predict(self, mutation_sets: Sequence[Tuple[str, ...]]) -> np.ndarray:
        return np.array([self.predict_one(ms) for ms in mutation_sets], dtype=float)


class NabuProteinModel:
    """
    Full NABU Protein Design Engine.
    Combines Component Additive Memory with Pairwise Residual Epistasis and Empirical Bayesian Shrinkage.
    """

    def __init__(
        self,
        shrinkage_prior: float = 1.0,
        aggregation: Literal["sum", "average"] = "sum"
    ):
        self.shrinkage_prior = shrinkage_prior
        self.aggregation = aggregation
        self.base_model = ComponentMemoryModel(shrinkage_prior=shrinkage_prior)
        self.pair_residuals: Dict[Tuple[str, str], Tuple[float, float]] = {}

    def fit(self, mutation_sets: Sequence[Tuple[str, ...]], fitness_scores: Sequence[float]) -> NabuProteinModel:
        scores = np.asarray(fitness_scores, dtype=float)
        
        # 1. Fit component additive memory
        self.base_model.fit(mutation_sets, scores)
        base_preds = self.base_model.predict(mutation_sets)
        residuals = scores - base_preds

        # 2. Extract and shrink pairwise residual interactions
        all_pairs: Set[Tuple[str, str]] = {
            p for ms in mutation_sets for p in combinations(ms, 2)
        }
        
        self.pair_residuals = {}
        for p in all_pairs:
            pair_set = set(p)
            indices = [i for i, ms in enumerate(mutation_sets) if pair_set.issubset(ms)]
            n = len(indices)
            if n > 0:
                pair_res_mean = float(np.mean(residuals[indices]))
                shrinkage = n / (n + self.shrinkage_prior)
                self.pair_residuals[p] = (pair_res_mean * shrinkage, shrinkage)

        return self

    def predict_one(self, mutation_set: Tuple[str, ...]) -> float:
        base_val = self.base_model.predict_one(mutation_set)
        
        pair_data = [
            self.pair_residuals[p]
            for p in combinations(mutation_set, 2)
            if p in self.pair_residuals
        ]
        
        if not pair_data:
            return base_val
            
        if self.aggregation == "sum":
            # Epistatic Potts / Ising energy model summation
            return base_val + sum(val * weight for val, weight in pair_data)
        else:
            # Weighted average consensus
            total_weighted_residual = sum(val * weight for val, weight in pair_data)
            total_weight = sum(weight for val, weight in pair_data)
            return base_val + (total_weighted_residual / total_weight if total_weight > 0 else 0.0)

    def predict(self, mutation_sets: Sequence[Tuple[str, ...]]) -> np.ndarray:
        return np.array([self.predict_one(ms) for ms in mutation_sets], dtype=float)


def evaluate_metrics(
    y_true: Sequence[float],
    y_pred: Sequence[float],
    candidate_budget: int = 50,
    top_quantile: float = 0.95
) -> Dict[str, float]:
    """
    Computes rigorous protein design benchmark metrics:
    - Spearman rank correlation
    - Top-1 percentile in ground truth
    - Top-5 mean percentile in ground truth
    - Enrichment factor of Top-50 in Top-5% ground truth
    - Normalized Top-1 regret
    """
    y_true_arr = np.asarray(y_true, dtype=float)
    y_pred_arr = np.asarray(y_pred, dtype=float)
    
    order = np.argsort(-y_pred_arr)
    sorted_y_true = y_true_arr[order]
    
    top1_true = sorted_y_true[0]
    top5_true = sorted_y_true[:5]
    top50_true = sorted_y_true[:candidate_budget]
    
    best_true = np.max(y_true_arr)
    worst_true = np.min(y_true_arr)
    cut = np.quantile(y_true_arr, top_quantile)
    
    spearman_corr = float(spearmanr(y_pred_arr, y_true_arr).statistic)
    top1_pct = percentile_rank(y_true_arr, top1_true)
    top5_mean_pct = float(np.mean([percentile_rank(y_true_arr, x) for x in top5_true]))
    
    baseline_prob = 1.0 - top_quantile
    enrichment = float(np.mean(top50_true >= cut) / baseline_prob)
    
    denom = best_true - worst_true
    regret = float((best_true - top1_true) / denom) if denom > 0 else 0.0
    
    return {
        "spearman": spearman_corr,
        "top1_percentile": top1_pct,
        "top5_mean_percentile": top5_mean_pct,
        "top50_enrichment": enrichment,
        "normalized_top1_regret": regret,
    }
