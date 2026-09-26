# NABU Phase-2D Root-Cause / Identifiability Audit

Date: 2026-09-26
Branch: exp/nabu-root-cause-audit
Dataset: revealed FLIP2 TrpB two-to-many only
NucB: MUST remain untouched
Phase 3: CLOSED

## Purpose

This experiment is diagnostic, not an architecture-selection benchmark.

The objective is to determine why low-order evidence (WT + singles + doubles) stops transferring to higher-order variants, and to separate:

1. missing higher-order information;
2. representation / estimation failure;
3. context-free interaction transfer failure;
4. composition / scaling failure.

No new reranker, blend weight, threshold search, or NucB evaluation is allowed in this audit.

## Frozen source contract

Official FLIP2 TrpB source:
- Zenodo record: 18433203
- file: trpb/two_to_many.csv.gz
- expected MD5: a611408d2db06907a023ddc8a1d94c12

Expected split:
- fit/train: 8,633
- validation: 2,158 (not used)
- test: 217,507

The test labels are intentionally revealed because TrpB is already a development dataset.

## Exact interaction decomposition

For a mutation set S, use a Möbius / epistasis decomposition.

WT term:
e0 = y(WT)

Single mutation:
e1(a) = y(a) - e0

Pair:
e2(a,b) = y(a,b) - e0 - e1(a) - e1(b)

Exact order-2 reconstruction:
O2(S) = e0 + sum e1 + sum e2

For a fully lower-order-supported triple:
e3(a,b,c) = y(a,b,c) - O2(a,b,c)

For a fully order-3-supported quad:
O3(a,b,c,d) = O2(a,b,c,d) + sum over its four triple e3 terms

and:
e4(a,b,c,d) = y(a,b,c,d) - O3(a,b,c,d)

These are measured diagnostic quantities, not deployable predictions.

## Arms

The audit must report:

- current NABU B2;
- current frozen B3 raw-pair score;
- exact singleton reconstruction O1;
- exact pair reconstruction O2;
- pair-identity shuffled O2 control;
- exact O3 diagnostic on supported quads.

## Coverage

Every oracle metric must include exact support counts.

No unsupported row may be silently imputed into an oracle arm.

B2/B3 must keep their existing fallback semantics.

## Required diagnostics

### A. Lower-order information transfer
For triples and quads separately:
- Spearman
- NDCG
- Top-1% hits
- RMSE
- MAE
- Pearson
- support fraction

### B. Interaction-order magnitude
Report:
- e3 mean/std/RMS/quantiles on supported triples;
- e4 mean/std/RMS/quantiles on supported quads;
- residual RMS divided by target standard deviation.

### C. Pair context stability
For every supported triple and each of its three pairs:
- base pair interaction = e2(pair)
- conditional pair interaction = e2(pair) + e3(triple)

Report:
- sign-flip rate;
- Spearman(base pair, conditional pair);
- Pearson(base pair, conditional pair).

### D. B2/B3 representation diagnostics
Report:
- shared mutation count between exact singleton effects and B2 main memory;
- correlation and scale ratio between B2 main effects and exact singleton effects;
- shared pair count between exact e2 and frozen B3 pair memory;
- correlation and scale ratio between B3 pair contribution and exact e2.

### E. Error attribution
On rows with O2 support, report correlation between:
- true O2 residual and B2 signed error;
- true O2 residual and B3 signed error.

For quads with O3 support, also report correlation between:
- exact e4 and B2 signed error;
- exact e4 and B3 signed error.

## Deterministic negative control

Pair effects are permuted across exact pair identities using seed 161 while preserving the empirical pair-effect distribution.

If exact pair identity is genuinely useful under context-free composition, exact O2 should outperform shuffled O2.

## Predeclared evidence flags

These flags are diagnostic labels, not scientific universal laws.

- MATERIAL_RANK_GAIN = absolute Spearman gain >= 0.10
- NEGLIGIBLE_PAIR_IDENTITY_GAIN = exact O2 Spearman - shuffled O2 Spearman <= 0.02
- MATERIAL_RESIDUAL = residual RMS / target SD >= 0.50
- MATERIAL_SIGN_INSTABILITY = conditional pair sign-flip rate >= 0.20
- REPRESENTATION_GAP = oracle Spearman - learned counterpart Spearman >= 0.10
- MATERIAL_B2_OR_B3_ERROR_ALIGNMENT = absolute Spearman(error, exact omitted residual) >= 0.50

## Root-cause decision tree

The script must not alter this logic after observing results.

1. If exact O2 is strong but B3 is worse by REPRESENTATION_GAP, flag:
   PAIR_REPRESENTATION_OR_ESTIMATION_FAILURE.

2. If O2 residual is MATERIAL_RESIDUAL and pair conditional interactions show MATERIAL_SIGN_INSTABILITY, flag:
   CONTEXT_DEPENDENT_HIGHER_ORDER_INTERACTION.

3. On quads, if exact O3 improves Spearman over exact O2 by MATERIAL_RANK_GAIN, flag:
   MISSING_THIRD_ORDER_INFORMATION_IS_CAUSAL_FOR_QUADS.

4. If quad e4 is MATERIAL_RESIDUAL after exact O3 reconstruction, flag:
   FOURTH_ORDER_INTERACTION_IS_MATERIAL.

5. If exact O2 does not materially beat shuffled O2, flag:
   CONTEXT_FREE_PAIR_IDENTITY_NOT_SUFFICIENT.

6. If B2 error strongly aligns with the exact omitted O2/O3 residual, flag:
   B2_FAILURE_TRACKS_MISSING_HIGHER_ORDER_TERMS.

The final JSON must list all triggered causes and select the primary cause by the following fixed priority:

1. MISSING_THIRD_ORDER_INFORMATION_IS_CAUSAL_FOR_QUADS
2. FOURTH_ORDER_INTERACTION_IS_MATERIAL
3. CONTEXT_DEPENDENT_HIGHER_ORDER_INTERACTION
4. PAIR_REPRESENTATION_OR_ESTIMATION_FAILURE
5. CONTEXT_FREE_PAIR_IDENTITY_NOT_SUFFICIENT
6. B2_FAILURE_TRACKS_MISSING_HIGHER_ORDER_TERMS
7. ROOT_CAUSE_NOT_RESOLVED_BY_THIS_AUDIT

Multiple causes may coexist and must not be hidden by the primary label.

## Interpretation boundary

A failure of exact context-free O2 does NOT prove that no model can ever infer higher-order fitness from low-order data.

It DOES show that a context-free WT + singles + pairs decomposition is insufficient on the tested region.

If measured triple terms materially recover quad ranking, those third-order coefficients are not directly identifiable from WT/single/double labels alone without an additional prior, representation, or additional measurements.

## Outputs

The run must save:

- ROOT_CAUSE_REPORT.json
- ROW_LEVEL_AUDIT.csv.gz
- INTERACTION_COMPONENTS.json
- RUN_MANIFEST.json

All outputs are development diagnostics and must be preserved even if the audit gives an inconvenient result.
