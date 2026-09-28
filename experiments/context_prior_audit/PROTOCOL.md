# NABU Frozen Context Incremental-Information Audit

Date: 2026-09-28
Branch: exp/nabu-root-cause-audit
Status: preregistered diagnostic; no architecture change

## Purpose

Stage 1-3 established that:
- implementation/data semantics are not the primary explanation;
- low-order identity representation is directionally recoverable;
- exact context-free pair terms are not portable across higher-order backgrounds;
- robust-asinh scale correction improves some rankings but does not remove context dependence.

This audit asks one narrower question before any new architecture is designed:

Does the already-frozen ESM contextual prior contain higher-order information that is independent of:
1. B2,
2. the pathological clean-additive baseline, and
3. mutation order?

No model is retrained in this audit.

## Immutable prediction inputs

### TrpB
- workflow run: 36243991521
- final artifact ID: 10907556059
- artifact name: nabu-v9-trpb-blind-result
- artifact digest: sha256:ea4fe1342c2d561d8574e21da1854aafa0549f2abebf132c1ada557867972a23
- prediction architecture was frozen before reveal.
- source MD5: a611408d2db06907a023ddc8a1d94c12

### IRED
- workflow run: 36237040369
- artifact ID: 10904493017
- artifact name: nabu-v9-contextual-esm-candidate
- artifact digest: sha256:a52c3338b8eaa01abe463484cb8bd376eaf3a7d2c70b80b0ab05919c62421fcc
- source SHA256 recorded by the artifact:
  aa45a2f85fb1af87b6b0e86397b3f8f292a061dc574be2536673b5bd490b0e74

NucB must not be downloaded or read.
Phase 3 remains closed.

## Primary analysis rows

Use only rows where:
complete_clean_main == True

Reason:
the contextual readout is only an active contextual correction on those rows.
Fallback rows where Context == B2 are excluded from the primary incremental-information test.

Minimum active rows required per dataset: 500.

## Variables

For every active row define:

T = revealed target
B = B2_MAIN_ONLY
A = clean additive baseline built from WT + measured singleton effects
P = frozen contextual residual prediction
C = A + P = frozen contextual final score

For IRED:
- A is the existing clean_additive_baseline_or_b2_fallback column on active rows.
- P is the existing contextual_residual column.
- C must equal V9_CONTEXTUAL_ESM_RESIDUAL within 1e-10.

For TrpB:
- reconstruct A from the pinned low-order fit source using exactly the historical rule:
  A = wt_target + sum(singleton_mean[mutation] - wt_target)
- P = V9_CONTEXTUAL_ESM_RESIDUAL - A
- C is the frozen V9_CONTEXTUAL_ESM_RESIDUAL.

The TrpB reconstruction must reproduce the historical complete_clean_main mask exactly.
If it does not, stop the audit.

## Primary conditional-information statistic

Use rank residualization.

For each dataset:
1. rank T, P, A, B with average ranks;
2. build nuisance matrix:
   intercept
   + rank(A)
   + rank(B)
   + one-hot mutation order (drop first category);
3. regress rank(T) on the nuisance matrix by ordinary least squares;
4. regress rank(P) on the same nuisance matrix;
5. compute Pearson correlation between the two residual vectors.

This is the primary conditional rank association:
rho_context_given_B2_additive_order

It asks whether frozen Context contains ordering information not already explained by B2, clean-additive cancellation, or mutation order.

## Permutation control

Use 1,000 deterministic permutations.

For each permutation:
- shuffle P only within mutation-order strata;
- seed stream derived from master seed 161;
- repeat the same rank-residualization procedure.

Report:
- real conditional association;
- permutation mean/std;
- one-sided empirical p-value:
  (1 + count(permuted >= real)) / 1001
- real minus permutation mean.

## Secondary direct-correction statistic

Define B2 error:
E_B = T - B

Compute Spearman(P, E_B), after residualizing rank(P) and rank(E_B)
against rank(A) + mutation-order one-hot.

This is secondary evidence only.

## Existing readout effectiveness

On exactly the same active rows report for B and C:
- Spearman
- NDCG
- Top-1% hits
- normalized regret

This is not used to decide whether P contains information.
It is used to distinguish "signal exists" from "current composition uses it correctly."

## Predeclared decisions

### INCREMENTAL_CONTEXT_SIGNAL_SUPPORTED
Only if BOTH TrpB and IRED satisfy:
- active rows >= 500
- primary conditional association >= 0.05
- empirical p <= 0.01
- real minus permutation mean >= 0.05

### FROZEN_CONTEXT_SIGNAL_PRESENT_BUT_READOUT_NOT_TRANSFERABLE
If INCREMENTAL_CONTEXT_SIGNAL_SUPPORTED is true AND
Context final-score Spearman is worse than B2 on at least one dataset's active rows.

### FROZEN_CONTEXT_SIGNAL_AND_READOUT_TRANSFER
If INCREMENTAL_CONTEXT_SIGNAL_SUPPORTED is true AND
Context final-score Spearman >= B2 on BOTH datasets' active rows.

### FROZEN_CONTEXT_NOT_CROSS_DATASET_INCREMENTAL
If the primary support gate fails on either dataset.

No result from this audit is allowed to imply that every possible protein language model prior succeeds or fails.

## Root-cause boundary

Even if Context has incremental information, do NOT declare the overall NABU root cause solved.

A positive result would establish:
- a frozen external contextual prior contains independent higher-order information;
- therefore the next architectural question is how to represent/use background-conditioned effects without destroying the B2 backbone.

A negative result would establish only:
- this frozen ESM residual prior is not reliable cross-dataset incremental evidence.

## Outputs

- CONTEXT_INCREMENTAL_MATRIX.json
- TRPB_CONTEXT_INCREMENTAL.json
- IRED_CONTEXT_INCREMENTAL.json
- RUN_MANIFEST.json
- OUTPUT_HASHES.json
