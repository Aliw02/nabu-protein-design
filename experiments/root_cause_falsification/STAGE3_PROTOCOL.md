# NABU Root-Cause Falsification — Stage 3 Scale-vs-Context Discrimination

Date: 2026-09-28
Branch: exp/nabu-root-cause-audit

## Motivation

Stage 2 supported three compatible explanations on both revealed development landscapes:
- H3 context-free pair non-portability
- H4 context-dependent higher-order interaction
- H5 observation-scale/composition nonlinearity

These are not yet independent causes. Stage 3 asks whether H5 explains most of H3/H4, or whether context dependence remains material after the prespecified scale correction.

Datasets:
- TrpB revealed development split
- GB1 revealed dense four-site landscape

NucB remains untouched.
Phase 3 remains closed.
No architecture is selected.

## Prespecified scale

Use robust-asinh only for the primary discrimination because it was the only prespecified Stage-2 transform that met the material cross-dataset criterion on matched order-3 rows.

Fit parameters use low-order targets only:
z = asinh((y - median_low_order) / MAD_low_order)

MAD fallback to low-order standard deviation only if MAD == 0.

Identity scale is the reference control.

Yeo-Johnson may be reported as secondary descriptive output but cannot determine the Stage-3 decision.

## Diagnostic A — Pair portability after scale correction

For identity and robust-asinh:
- derive WT, exact e1, exact e2 from low-order fit rows on that scale;
- evaluate exact O1, exact O2, and seed-161 shuffled-pair O2 on supported triples and quads.

For each dataset/order report:
- O1 Spearman
- O2 Spearman
- shuffled O2 Spearman
- O2-O1 gain
- O2-shuffled gain

Primary scale rescue criterion:
robust-asinh must improve O2 Spearman over identity by >= 0.10 on BOTH datasets for the same mutation order.

## Diagnostic B — Context dependence after scale correction

On exact-O2-supported triples, for identity and robust-asinh:
e3 = transformed_true_target - transformed_O2
conditional_pair = transformed_e2(pair) + e3(triple)

Report:
- sign-flip rate
- Spearman(base e2, conditional pair)
- Pearson(base e2, conditional pair)
- e3 RMS / transformed target SD
- Spearman(e3, O2)
- Pearson(e3, O2)

Define relative reduction from identity:
- sign_flip_reduction = 1 - robust_sign_flip / identity_sign_flip
- e3_ratio_reduction = 1 - robust_e3_ratio / identity_e3_ratio

Scale-mediated context criterion:
- sign-flip reduction >= 50%
- e3 RMS/target-SD reduction >= 50%
- BOTH on TrpB and GB1

Persistent context criterion:
after robust-asinh:
- sign-flip rate >= 0.20
- e3 RMS/target SD >= 0.50
- BOTH on TrpB and GB1

Possible outcomes:
1. SCALE_MEDIATES_MOST_CONTEXT_DEPENDENCE
2. CONTEXT_DEPENDENCE_PERSISTS_AFTER_SCALE_CORRECTION
3. MIXED_OR_UNRESOLVED

If both scale-mediated and persistent criteria happen to hold, classify MIXED_OR_UNRESOLVED.

## Diagnostic C — Cross-order cancellation transfer

Purpose: test whether one simple cancellation law learned from triples transfers to quads.

On each dataset:
- use only exact-O2-supported triples to fit:
  residual3 = a + b * O2
  by ordinary least squares;
- freeze a,b;
- apply to exact-O2-supported quads:
  predicted_quad = O2_quad + a + b * O2_quad

Controls:
- raw O2 quad prediction
- B2 quad prediction
- seed-161 shuffled triple residual pairing, fit same affine rule, then apply to quads

Report on quads:
- Spearman
- Pearson
- RMSE
- NDCG where defined

Cross-order cancellation transfer is SUPPORTED only if on BOTH datasets:
- affine-triple-law quad Spearman >= raw O2 quad Spearman + 0.10
- affine-triple-law quad Spearman >= shuffled-law quad Spearman + 0.10

This remains a revealed-data diagnostic, not a deployable model claim.

## Diagnostic D — Measured third-order contribution to quads

For quads whose four constituent triples are measured and exact-O2 supported:
- derive measured e3 for each constituent triple;
- O3 = O2 + sum constituent e3;
- e4 = true quad - O3

Run on identity and robust-asinh.

Report:
- supported quad count
- O2 Spearman
- O3 Spearman
- O3-O2 gain
- e4 RMS / target SD
- Spearman(e4, O3)

Interpretation:
- if robust-asinh makes O3 unnecessary and e4 small, scale explains most apparent higher-order correction;
- if large e4 remains after robust-asinh, scale alone is insufficient.

Material remaining fourth-order residual:
e4 RMS / target SD >= 0.50.

## Decision language

Stage 3 may discriminate SCALE vs PERSISTENT CONTEXT, but it must not declare the entire NABU root cause unless all remaining competing explanations are eliminated.

Global output remains ROOT_CAUSE_NOT_YET_IDENTIFIED unless one mechanism explains the others under the preregistered discrimination tests.

## Outputs

- STAGE3_MATRIX.json
- SCALE_CONTEXT_DIAGNOSTICS.json
- CROSS_ORDER_TRANSFER.json
- QUAD_DECOMPOSITION.json
- RUN_MANIFEST.json
- OUTPUT_HASHES.json
