#!/usr/bin/env python3
"""Lightning Protein/DMS POC V1 ablation.

This script reuses only the already-frozen mutant->fold assignment from POC V1.
For each fold the predictor receives:
  * visible rows with DMS_score/DMS_score_bin;
  * hidden descriptors (mutant, position, WT, ALT, context) WITHOUT hidden truth.
Predictions are serialized and SHA-256 frozen before hidden truth is merged.

No amino-acid property tables, domain vocabulary, embeddings, or hidden-label
threshold tuning are used.
"""
from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
import os
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from scipy.stats import spearmanr, wilcoxon
from sklearn.metrics import average_precision_score, roc_auc_score

ROOT = Path(__file__).resolve().parent
DATASET = Path('/mnt/data/kcnh2_subset.csv')
PRIOR_PRED = Path('/mnt/data/lightning_protein_poc_v1/predictions_revealed.csv')
SEED = 161
N_BOOT = 5000

AA = set('ACDEFGHIKLMNPQRSTVWY')


def parse_mutant(mutant: str) -> Tuple[str, int, str]:
    return mutant[0], int(mutant[1:-1]), mutant[-1]


def read_csv(path: Path) -> List[dict]:
    with path.open(newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: List[dict], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def context_for(pos: int, wt_by_pos: Dict[int, str]) -> str:
    left = wt_by_pos.get(pos - 1, '^')
    center = wt_by_pos[pos]
    right = wt_by_pos.get(pos + 1, '$')
    return f'{left}[{center}]{right}'


def mean(values):
    values = list(values)
    return float(sum(values) / len(values)) if values else float('nan')


def build_relative_resonance(visible: List[dict], positions: List[int]):
    """Return mutual row-relative resonance learned from visible response profiles.

    Raw affinity is Pearson correlation between position response profiles across
    amino-acid alternatives shared in the visible set (>=3 shared alternatives).
    Relative attention follows the same row-relative form used in the current
    Lightning A5 work: exp((corr-row_max)/row_positive_std), followed by mutual
    geometric agreement. Negative/non-finite correlations do not support Dream.
    """
    by_pos_alt = defaultdict(dict)
    for r in visible:
        by_pos_alt[r['pos']][r['alt']] = r['score']

    n = len(positions)
    corr = np.full((n, n), np.nan, dtype=float)
    for i, p in enumerate(positions):
        corr[i, i] = 1.0
        for j in range(i + 1, n):
            q = positions[j]
            shared = sorted(set(by_pos_alt[p]) & set(by_pos_alt[q]))
            if len(shared) < 3:
                continue
            x = np.array([by_pos_alt[p][a] for a in shared], dtype=float)
            y = np.array([by_pos_alt[q][a] for a in shared], dtype=float)
            if np.std(x) <= 1e-12 or np.std(y) <= 1e-12:
                c = 0.0
            else:
                c = float(np.corrcoef(x, y)[0, 1])
                if not math.isfinite(c):
                    c = 0.0
            corr[i, j] = corr[j, i] = c

    relative = np.zeros((n, n), dtype=float)
    for i in range(n):
        row = corr[i].copy()
        row[i] = -np.inf
        positive = row[np.isfinite(row) & (row > 0.0)]
        if positive.size == 0:
            continue
        maximum = float(np.max(positive))
        scale = float(np.std(positive))
        if (not math.isfinite(scale)) or scale <= 1e-12:
            scale = max(abs(maximum), 1.0)
        valid = np.isfinite(row) & (row > 0.0)
        exponent = np.maximum((row[valid] - maximum) / scale, -700.0)
        relative[i, valid] = np.exp(exponent)

    mutual = np.sqrt(relative * relative.T)
    np.fill_diagonal(mutual, 1.0)
    return corr, relative, mutual


def predict_fold(fold: int, visible: List[dict], hidden_desc: List[dict], positions: List[int], wt_by_pos: Dict[int, str]):
    # Strict API boundary: this function has no hidden score/bin parameter.

    # Visible-only statistics.
    global_mean = mean(r['score'] for r in visible)
    pos_scores = defaultdict(list)
    alt_scores = defaultdict(list)
    pos_alt_score = {}
    for r in visible:
        pos_scores[r['pos']].append(r['score'])
        alt_scores[r['alt']].append(r['score'])
        pos_alt_score[(r['pos'], r['alt'])] = r['score']
    pos_mean = {p: mean(v) for p, v in pos_scores.items()}
    alt_mean = {a: mean(v) for a, v in alt_scores.items()}

    raw_corr, relative, mutual = build_relative_resonance(visible, positions)
    pindex = {p: i for i, p in enumerate(positions)}

    preds = []
    for h in hidden_desc:
        p = h['pos']
        alt = h['alt']

        b0 = global_mean
        b1 = pos_mean.get(p, global_mean)
        alt_global = alt_mean.get(alt, global_mean)

        # Position-centered ALT memory: learn the ALT effect as the mean residual
        # across other visible positions. This is the exact POC V1 CORE_ADDITIVE
        # semantics and avoids confounding ALT effects with position baselines.
        alt_residuals = []
        for q in positions:
            s = pos_alt_score.get((q, alt))
            if s is None or q not in pos_mean:
                continue
            alt_residuals.append(float(s - pos_mean[q]))
        mean_alt_residual = mean(alt_residuals) if alt_residuals else 0.0
        b2 = b1 + mean_alt_residual

        # Relative-resonance memory transfer from the same ALT at other visible positions.
        supports = []
        i = pindex[p]
        for q in positions:
            if q == p:
                continue
            j = pindex[q]
            w = float(mutual[i, j])
            if w <= 0.0 or not math.isfinite(w):
                continue
            s = pos_alt_score.get((q, alt))
            if s is None:
                continue
            residual = float(s - pos_mean[q])
            supports.append((q, w, residual))

        if supports:
            denom = sum(w for _, w, _ in supports)
            transfer = sum(w * residual for _, w, residual in supports) / denom
            b3 = b1 + transfer
            status = 'WARM'
        else:
            # No positive resonance support: preserve Memory estimate; Dream state is DORMANT.
            b3 = b2
            status = 'DORMANT'

        typed = {
            'CORE': f'position:{p}',
            'IDENTITY_ANCHOR': f'{h["wt"]}{p}',
            'QUALIFIER': f'alt:{alt}',
            'QUANTITY': 'mutation_count:1',
            'CONTEXT': h['context'],
        }
        preds.append({
            **h,
            'B0_GLOBAL_MEAN': b0,
            'B1_A4_CORE_POSITION': b1,
            'B1_ALT_GLOBAL_REFERENCE': alt_global,
            'B2_A4_PLUS_MEMORY': b2,
            'B3_PLUS_RELATIVE_RESONANCE': b3,
            'B4_DREAM_WARM_ONLY': b3 if status == 'WARM' else None,
            'B5_FULL_DREAM_DORMANT': b3,
            'dream_status': status,
            'supporting_positions': [q for q, _, _ in supports],
            'positive_resonance_support_count': len(supports),
            'max_mutual_resonance': max((w for _, w, _ in supports), default=0.0),
            'typed_contract': typed,
        })

    # Freeze predictions BEFORE truth reveal.
    freeze_payload = {
        'version': 'LIGHTNING_PROTEIN_POC_ABLATION_V1',
        'fold': fold,
        'hidden_count': len(preds),
        'predictor_input_contract': {
            'visible_has_truth': True,
            'hidden_descriptor_fields': ['mutant', 'fold', 'pos', 'wt', 'alt', 'context'],
            'hidden_truth_fields': [],
        },
        'predictions': preds,
    }
    encoded = json.dumps(freeze_payload, sort_keys=True, separators=(',', ':'), ensure_ascii=False).encode('utf-8')
    sha = hashlib.sha256(encoded).hexdigest()
    freeze_path = ROOT / 'freezes' / f'fold_{fold}_pre_reveal.json'
    freeze_path.parent.mkdir(parents=True, exist_ok=True)
    freeze_path.write_bytes(json.dumps(freeze_payload, ensure_ascii=False, indent=2, sort_keys=True).encode('utf-8'))

    for p in preds:
        p['freeze_sha256'] = sha
    return preds, sha


def metric_vector(rows, method):
    selected = [(float(r[method]), float(r['score']), int(r['bin'])) for r in rows if r.get(method) is not None]
    if not selected:
        return np.array([]), np.array([]), np.array([])
    pred = np.array([x[0] for x in selected], dtype=float)
    truth = np.array([x[1] for x in selected], dtype=float)
    binary = np.array([x[2] for x in selected], dtype=int)
    return pred, truth, binary


def ranking_group_metrics(rows, method):
    groups = defaultdict(list)
    for r in rows:
        if r.get(method) is None:
            continue
        groups[(r['fold'], r['pos'])].append(r)

    per_group = []
    for key, g in sorted(groups.items()):
        g_pred = sorted(g, key=lambda r: (-float(r[method]), r['mutant']))
        g_truth = sorted(g, key=lambda r: (-float(r['score']), r['mutant']))
        best_mutant = g_truth[0]['mutant']
        ranks = {r['mutant']: i + 1 for i, r in enumerate(g_pred)}
        best_rank = ranks[best_mutant]
        top1 = g_pred[:1]
        top3 = g_pred[:3]
        top5 = g_pred[:5]
        per_group.append({
            'fold': key[0], 'pos': key[1], 'candidate_count': len(g),
            'top1_score': mean(r['score'] for r in top1),
            'top3_score': mean(r['score'] for r in top3),
            'top5_score': mean(r['score'] for r in top5),
            'top1_positive_rate': mean(r['bin'] for r in top1),
            'top3_positive_rate': mean(r['bin'] for r in top3),
            'top5_positive_rate': mean(r['bin'] for r in top5),
            'best_hidden_rank': best_rank,
            'hit_best_at1': int(best_rank <= 1),
            'hit_best_at3': int(best_rank <= 3),
            'hit_best_at5': int(best_rank <= 5),
            'top1_regret': float(g_truth[0]['score'] - top1[0]['score']),
        })
    return per_group


def method_metrics(rows, method):
    pred, truth, binary = metric_vector(rows, method)
    coverage = len(pred) / len(rows)
    rho = float(spearmanr(pred, truth).statistic) if len(pred) >= 3 else float('nan')
    auc = float(roc_auc_score(binary, pred)) if len(set(binary.tolist())) == 2 else float('nan')
    ap = float(average_precision_score(binary, pred)) if len(pred) else float('nan')
    groups = ranking_group_metrics(rows, method)
    base_rate = mean(r['bin'] for r in rows if r.get(method) is not None)
    top1_pos = mean(g['top1_positive_rate'] for g in groups)
    top3_pos = mean(g['top3_positive_rate'] for g in groups)
    top5_pos = mean(g['top5_positive_rate'] for g in groups)
    return {
        'method': method,
        'coverage': coverage,
        'n_predictions': len(pred),
        'spearman': rho,
        'roc_auc_bin': auc,
        'avg_precision_bin': ap,
        'base_positive_rate': base_rate,
        'mean_top1_hidden_score': mean(g['top1_score'] for g in groups),
        'mean_top3_hidden_score': mean(g['top3_score'] for g in groups),
        'mean_top5_hidden_score': mean(g['top5_score'] for g in groups),
        'top1_positive_rate': top1_pos,
        'top3_positive_rate': top3_pos,
        'top5_positive_rate': top5_pos,
        'top1_enrichment_vs_base': top1_pos / base_rate if base_rate else float('nan'),
        'top3_enrichment_vs_base': top3_pos / base_rate if base_rate else float('nan'),
        'top5_enrichment_vs_base': top5_pos / base_rate if base_rate else float('nan'),
        'median_best_hidden_rank': float(np.median([g['best_hidden_rank'] for g in groups])),
        'mean_best_hidden_rank': mean(g['best_hidden_rank'] for g in groups),
        'hit_best_at1': mean(g['hit_best_at1'] for g in groups),
        'hit_best_at3': mean(g['hit_best_at3'] for g in groups),
        'hit_best_at5': mean(g['hit_best_at5'] for g in groups),
        'mean_top1_regret': mean(g['top1_regret'] for g in groups),
        'group_count': len(groups),
    }, groups


def per_position_spearman(rows, method):
    output = {}
    by_pos = defaultdict(list)
    for r in rows:
        if r.get(method) is not None:
            by_pos[r['pos']].append(r)
    for p, g in sorted(by_pos.items()):
        pred = [r[method] for r in g]
        truth = [r['score'] for r in g]
        output[p] = float(spearmanr(pred, truth).statistic) if len(g) >= 3 else float('nan')
    return output


def exact_sign_flip_pvalue(differences):
    d = np.array([x for x in differences if math.isfinite(x) and abs(x) > 1e-15], dtype=float)
    if d.size == 0:
        return 1.0
    observed = abs(float(np.mean(d)))
    if d.size <= 20:
        vals = []
        for signs in itertools.product((-1.0, 1.0), repeat=d.size):
            vals.append(abs(float(np.mean(d * np.array(signs)))))
        vals = np.array(vals)
        return float((np.sum(vals >= observed - 1e-15)) / len(vals))
    rng = np.random.default_rng(SEED)
    count = 0
    n = 20000
    for _ in range(n):
        signs = rng.choice([-1.0, 1.0], size=d.size)
        if abs(float(np.mean(d * signs))) >= observed - 1e-15:
            count += 1
    return (count + 1) / (n + 1)


def cluster_bootstrap_delta(rows, method_a, method_b, metric='spearman', n_boot=N_BOOT):
    positions = sorted({r['pos'] for r in rows})
    by_pos = {p: [r for r in rows if r['pos'] == p] for p in positions}
    rng = np.random.default_rng(SEED)
    deltas = []
    for _ in range(n_boot):
        sampled = rng.choice(positions, size=len(positions), replace=True)
        boot = []
        # Duplicate clusters are intentionally duplicated to weight the resample.
        for draw_idx, p in enumerate(sampled):
            for r in by_pos[int(p)]:
                rr = dict(r)
                rr['_bootstrap_draw'] = draw_idx
                boot.append(rr)
        pa, ta, ba = metric_vector(boot, method_a)
        pb, tb, bb = metric_vector(boot, method_b)
        if metric == 'spearman':
            va = float(spearmanr(pa, ta).statistic)
            vb = float(spearmanr(pb, tb).statistic)
        elif metric == 'auc':
            if len(set(ba.tolist())) < 2 or len(set(bb.tolist())) < 2:
                continue
            va = float(roc_auc_score(ba, pa))
            vb = float(roc_auc_score(bb, pb))
        else:
            raise ValueError(metric)
        if math.isfinite(va) and math.isfinite(vb):
            deltas.append(vb - va)
    arr = np.array(deltas, dtype=float)
    return {
        'n_boot': int(arr.size),
        'delta_mean': float(np.mean(arr)),
        'ci95_low': float(np.quantile(arr, 0.025)),
        'ci95_high': float(np.quantile(arr, 0.975)),
        'fraction_delta_gt_0': float(np.mean(arr > 0)),
    }


def main():
    raw = read_csv(DATASET)
    prior = read_csv(PRIOR_PRED)
    fold_map = {r['mutant']: int(r['fold']) for r in prior}

    rows = []
    wt_by_pos = {}
    for r in raw:
        wt, pos, alt = parse_mutant(r['mutant'])
        wt_by_pos[pos] = wt
        rows.append({
            'mutant': r['mutant'], 'wt': wt, 'pos': pos, 'alt': alt,
            'score': float(r['DMS_score']), 'bin': int(r['DMS_score_bin']),
            'fold': fold_map[r['mutant']],
        })
    positions = sorted(wt_by_pos)
    for r in rows:
        r['context'] = context_for(r['pos'], wt_by_pos)

    # Produce pre-reveal predictions and freeze each fold before truth merge.
    predicted = []
    freezes = []
    for fold in sorted({r['fold'] for r in rows}):
        visible_fold = [dict(r) for r in rows if r['fold'] != fold]
        hidden_desc_fold = [
            {k: r[k] for k in ('mutant', 'fold', 'pos', 'wt', 'alt', 'context')}
            for r in rows if r['fold'] == fold
        ]
        p, sha = predict_fold(fold, visible_fold, hidden_desc_fold, positions, wt_by_pos)
        predicted.extend(p)
        freezes.append({'fold': fold, 'hidden_count': len(p), 'sha256': sha})

    truth = {r['mutant']: r for r in rows}
    revealed = []
    for p in predicted:
        t = truth[p['mutant']]
        rr = dict(p)
        rr['DMS_score'] = t['score']
        rr['DMS_score_bin'] = t['bin']
        # Internal standardized names for metrics.
        rr['score'] = t['score']
        rr['bin'] = t['bin']
        revealed.append(rr)

    methods = [
        'B0_GLOBAL_MEAN',
        'B1_A4_CORE_POSITION',
        'B1_ALT_GLOBAL_REFERENCE',
        'B2_A4_PLUS_MEMORY',
        'B3_PLUS_RELATIVE_RESONANCE',
        'B4_DREAM_WARM_ONLY',
        'B5_FULL_DREAM_DORMANT',
    ]
    metrics = []
    group_rows = []
    groups_by_method = {}
    for m in methods:
        mm, gg = method_metrics(revealed, m)
        metrics.append(mm)
        groups_by_method[m] = gg
        for g in gg:
            group_rows.append({'method': m, **g})

    # Component delta table.
    metric_by_name = {m['method']: m for m in metrics}
    chain = [
        ('CORE_CONTEXT', 'B0_GLOBAL_MEAN', 'B1_A4_CORE_POSITION'),
        ('MEMORY', 'B1_A4_CORE_POSITION', 'B2_A4_PLUS_MEMORY'),
        ('RELATIVE_RESONANCE', 'B2_A4_PLUS_MEMORY', 'B3_PLUS_RELATIVE_RESONANCE'),
        ('DREAM_WARM_FILTER', 'B3_PLUS_RELATIVE_RESONANCE', 'B4_DREAM_WARM_ONLY'),
        ('DORMANT_RETENTION', 'B4_DREAM_WARM_ONLY', 'B5_FULL_DREAM_DORMANT'),
    ]
    component_deltas = []
    for component, before, after in chain:
        a, b = metric_by_name[before], metric_by_name[after]
        component_deltas.append({
            'component': component, 'before': before, 'after': after,
            'delta_coverage': b['coverage'] - a['coverage'],
            'delta_spearman': b['spearman'] - a['spearman'],
            'delta_auc': b['roc_auc_bin'] - a['roc_auc_bin'],
            'delta_ap': b['avg_precision_bin'] - a['avg_precision_bin'],
            'delta_top1_score': b['mean_top1_hidden_score'] - a['mean_top1_hidden_score'],
            'delta_top3_score': b['mean_top3_hidden_score'] - a['mean_top3_hidden_score'],
            'delta_top5_score': b['mean_top5_hidden_score'] - a['mean_top5_hidden_score'],
            'delta_hit_best_at3': b['hit_best_at3'] - a['hit_best_at3'],
            'delta_top1_regret': b['mean_top1_regret'] - a['mean_top1_regret'],
        })

    # Exact reproduction audit against POC V1 predictions.
    prior_by_mut = {r['mutant']: r for r in prior}
    reproduction = {}
    comparisons = {
        'B1_A4_CORE_POSITION_vs_prior_B0_POSITION_MEAN': ('B1_A4_CORE_POSITION', 'B0_POSITION_MEAN'),
        'B1_ALT_GLOBAL_REFERENCE_vs_prior_B1_ALT_GLOBAL': ('B1_ALT_GLOBAL_REFERENCE', 'B1_ALT_GLOBAL'),
        'B2_A4_PLUS_MEMORY_vs_prior_CORE_ADDITIVE': ('B2_A4_PLUS_MEMORY', 'CORE_ADDITIVE'),
        'B5_FULL_DREAM_DORMANT_vs_prior_LIGHTNING_DREAM': ('B5_FULL_DREAM_DORMANT', 'LIGHTNING_DREAM'),
        'max_mutual_resonance_vs_prior': ('max_mutual_resonance', 'max_mutual_resonance'),
    }
    for label, (newk, oldk) in comparisons.items():
        diffs = []
        for r in revealed:
            diffs.append(abs(float(r[newk]) - float(prior_by_mut[r['mutant']][oldk])))
        reproduction[label] = {'max_abs_delta': max(diffs), 'mean_abs_delta': mean(diffs)}

    # Per-position paired analysis for the critical incremental gain: Memory -> Resonance/Full.
    pos_b2 = per_position_spearman(revealed, 'B2_A4_PLUS_MEMORY')
    pos_b5 = per_position_spearman(revealed, 'B5_FULL_DREAM_DORMANT')
    pos_rows = []
    diffs = []
    for p in positions:
        d = pos_b5[p] - pos_b2[p]
        diffs.append(d)
        pos_rows.append({
            'pos': p, 'wt': wt_by_pos[p],
            'B2_memory_spearman': pos_b2[p],
            'B5_full_spearman': pos_b5[p],
            'delta': d,
            'n_variants': sum(1 for r in revealed if r['pos'] == p),
        })
    sign_flip_p = exact_sign_flip_pvalue(diffs)

    # Top-3 paired group test across fold x position.
    g2 = {(g['fold'], g['pos']): g for g in groups_by_method['B2_A4_PLUS_MEMORY']}
    g5 = {(g['fold'], g['pos']): g for g in groups_by_method['B5_FULL_DREAM_DORMANT']}
    common = sorted(set(g2) & set(g5))
    top3_b2 = np.array([g2[k]['top3_score'] for k in common], dtype=float)
    top3_b5 = np.array([g5[k]['top3_score'] for k in common], dtype=float)
    if np.allclose(top3_b2, top3_b5):
        top3_w = {'statistic': 0.0, 'pvalue': 1.0}
    else:
        wr = wilcoxon(top3_b5, top3_b2, alternative='two-sided', zero_method='wilcox')
        top3_w = {'statistic': float(wr.statistic), 'pvalue': float(wr.pvalue)}

    bootstrap_spear = cluster_bootstrap_delta(revealed, 'B2_A4_PLUS_MEMORY', 'B5_FULL_DREAM_DORMANT', 'spearman')
    bootstrap_auc = cluster_bootstrap_delta(revealed, 'B2_A4_PLUS_MEMORY', 'B5_FULL_DREAM_DORMANT', 'auc')

    # Candidate-state audit.
    warm = sum(1 for r in revealed if r['dream_status'] == 'WARM')
    dormant = sum(1 for r in revealed if r['dream_status'] == 'DORMANT')
    dormant_rows = [r for r in revealed if r['dream_status'] == 'DORMANT']

    # Save tables.
    metric_fields = list(metrics[0].keys())
    write_csv(ROOT / 'ablation_metrics.csv', metrics, metric_fields)
    write_csv(ROOT / 'component_deltas.csv', component_deltas, list(component_deltas[0].keys()))
    write_csv(ROOT / 'per_position_spearman.csv', pos_rows, list(pos_rows[0].keys()))
    write_csv(ROOT / 'group_ranking_metrics.csv', group_rows, list(group_rows[0].keys()))

    pred_fields = [
        'mutant','fold','pos','wt','alt','context',
        'B0_GLOBAL_MEAN','B1_A4_CORE_POSITION','B1_ALT_GLOBAL_REFERENCE',
        'B2_A4_PLUS_MEMORY','B3_PLUS_RELATIVE_RESONANCE',
        'B4_DREAM_WARM_ONLY','B5_FULL_DREAM_DORMANT',
        'dream_status','supporting_positions','positive_resonance_support_count',
        'max_mutual_resonance','typed_contract','freeze_sha256','DMS_score','DMS_score_bin'
    ]
    out_rows = []
    for r in revealed:
        out = {k: r.get(k) for k in pred_fields}
        out['supporting_positions'] = json.dumps(out['supporting_positions'])
        out['typed_contract'] = json.dumps(out['typed_contract'], sort_keys=True)
        out_rows.append(out)
    write_csv(ROOT / 'predictions_revealed.csv', out_rows, pred_fields)

    significance = {
        'critical_comparison': 'B2_A4_PLUS_MEMORY -> B5_FULL_DREAM_DORMANT',
        'per_position_spearman_exact_sign_flip': {
            'positions': len(diffs),
            'mean_delta': mean(diffs),
            'positive_positions': sum(1 for d in diffs if d > 0),
            'negative_positions': sum(1 for d in diffs if d < 0),
            'zero_positions': sum(1 for d in diffs if abs(d) <= 1e-15),
            'two_sided_pvalue': sign_flip_p,
        },
        'cluster_bootstrap_spearman_delta_B5_minus_B2': bootstrap_spear,
        'cluster_bootstrap_auc_delta_B5_minus_B2': bootstrap_auc,
        'fold_position_top3_score_wilcoxon_B5_vs_B2': {
            'group_count': len(common), **top3_w,
            'mean_delta': float(np.mean(top3_b5 - top3_b2)),
            'median_delta': float(np.median(top3_b5 - top3_b2)),
        },
    }
    (ROOT / 'significance.json').write_text(json.dumps(significance, indent=2), encoding='utf-8')

    audit = {
        'version': 'LIGHTNING_PROTEIN_POC_ABLATION_V1',
        'dataset': DATASET.name,
        'rows': len(rows),
        'positions': positions,
        'same_frozen_fold_assignment_as_poc_v1': True,
        'seed': SEED,
        'predictor_hidden_input_fields': ['mutant','fold','pos','wt','alt','context'],
        'predictor_hidden_truth_fields': [],
        'hidden_truth_used_before_freeze': False,
        'amino_acid_property_table_used': False,
        'manual_protein_domain_vocabulary_used': False,
        'relative_resonance_source': 'visible-only Pearson response-profile correlations transformed by row-relative mutual attention',
        'freeze_records': freezes,
        'reproduction_against_poc_v1': reproduction,
        'dream_states': {'WARM': warm, 'DORMANT': dormant},
        'dormant_candidates': [
            {'mutant': r['mutant'], 'fold': r['fold'], 'DMS_score': r['score'], 'B5_score': r['B5_FULL_DREAM_DORMANT']}
            for r in dormant_rows
        ],
    }
    (ROOT / 'leakage_and_reproduction_audit.json').write_text(json.dumps(audit, indent=2), encoding='utf-8')

    summary = {
        'version': 'LIGHTNING_PROTEIN_POC_ABLATION_V1',
        'metrics': metrics,
        'component_deltas': component_deltas,
        'significance': significance,
        'audit': audit,
    }
    (ROOT / 'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')

    # Human-readable README.
    m = metric_by_name
    readme = f'''# Lightning Protein/DMS POC — Ablation V1\n\nDataset: 200 KCNH2 single-mutant DMS variants, positions 545–555.\n\n## Protocol\n\n- Reuses the exact frozen 5-fold mutant assignment from POC V1.\n- Hidden predictor input contains descriptors only; hidden DMS truth is absent.\n- Each fold prediction payload is SHA-256 frozen before truth merge.\n- No amino-acid property table, protein-domain vocabulary, embeddings, or hidden-label tuning.\n- Relative resonance is learned only from visible response profiles.\n\n## Arms\n\n- **B0_GLOBAL_MEAN** — visible global mean.\n- **B1_A4_CORE_POSITION** — typed CORE/context represented by the visible mean at the exact position.\n- **B2_A4_PLUS_MEMORY** — B1 plus position-centered ALT residual memory learned across visible positions.\n- **B3_PLUS_RELATIVE_RESONANCE** — B2 upgraded by mutual row-relative resonance transfer when visible same-ALT support exists; B2 fallback otherwise.\n- **B4_DREAM_WARM_ONLY** — same score as B3 but only WARM candidates are retained.\n- **B5_FULL_DREAM_DORMANT** — WARM candidates plus unsupported candidates retained DORMANT with the B2 estimate.\n\n## Critical result\n\nB2 Memory Spearman: **{m['B2_A4_PLUS_MEMORY']['spearman']:.6f}**  \nB5 Full Spearman: **{m['B5_FULL_DREAM_DORMANT']['spearman']:.6f}**  \nDelta: **{m['B5_FULL_DREAM_DORMANT']['spearman']-m['B2_A4_PLUS_MEMORY']['spearman']:+.6f}**\n\nB2 ROC-AUC: **{m['B2_A4_PLUS_MEMORY']['roc_auc_bin']:.6f}**  \nB5 ROC-AUC: **{m['B5_FULL_DREAM_DORMANT']['roc_auc_bin']:.6f}**  \nDelta: **{m['B5_FULL_DREAM_DORMANT']['roc_auc_bin']-m['B2_A4_PLUS_MEMORY']['roc_auc_bin']:+.6f}**\n\nB2 mean Top-3 hidden score: **{m['B2_A4_PLUS_MEMORY']['mean_top3_hidden_score']:.4f}**  \nB5 mean Top-3 hidden score: **{m['B5_FULL_DREAM_DORMANT']['mean_top3_hidden_score']:.4f}**\n\nPer-position Spearman sign-flip p-value (B5 vs B2): **{sign_flip_p:.6g}**.  \nCluster bootstrap 95% CI for Spearman delta: **[{bootstrap_spear['ci95_low']:.6f}, {bootstrap_spear['ci95_high']:.6f}]**.\n\nDream states: **{warm} WARM / {dormant} DORMANT**.\n\n## Interpretation boundary\n\nThis is an engineering ablation on one small DMS slice. It tests whether the Lightning-inspired representation/memory/resonance mechanics add blind predictive ranking signal. It does **not** establish biological validity, de-novo protein design performance, or laboratory efficacy.\n'''
    (ROOT / 'README.md').write_text(readme, encoding='utf-8')

    print(json.dumps({
        'metrics': metrics,
        'component_deltas': component_deltas,
        'significance': significance,
        'reproduction': reproduction,
        'dream_states': {'WARM': warm, 'DORMANT': dormant},
    }, indent=2))


if __name__ == '__main__':
    main()
