# NABU Phase-2D Root-Cause Falsification Protocol

Date: 2026-09-28
Branch: exp/nabu-root-cause-audit
Scope: revealed development datasets only (TrpB + IRED)
Untouched benchmark: NucB MUST NOT be downloaded, read, or evaluated.
Phase 3: CLOSED.

## Purpose

Stop the patch -> fail -> reinterpret loop.

This protocol does not select a new architecture. It tests competing explanations and records each as SUPPORTED, REJECTED, or UNRESOLVED.

No hypothesis becomes "the root cause" from one dataset or one metric.

## Stage 0 — Implementation / data invariants

Run independently on TrpB and IRED.

Required invariants:
1. The derived reference has one fixed length and valid amino-acid alphabet.
2. mutation_set(reference, reference) is empty.
3. Every fit sequence reconstructs exactly from reference + parsed mutation tokens.
4. Exact O1 reconstructs every WT and singleton fit target to absolute error <= 1e-12.
5. Exact O2 reconstructs every clean double fit target to absolute error <= 1e-12.
6. Exact O2 support bookkeeping equals the direct combinatorial support check.
7. Metric self-tests:
   - Spearman(y, y) == 1 within 1e-12.
   - Pearson(y, y) == 1 within 1e-12.
   - RMSE(y, y) == 0 within 1e-12.
8. Split counts and pinned source hashes must match the historical contracts.

If any required invariant fails:
- IMPLEMENTATION_OR_DATA_SEMANTICS = SUPPORTED
- stop biological/mechanistic interpretation for that dataset.

If all required invariants pass on both datasets:
- IMPLEMENTATION_OR_DATA_SEMANTICS = REJECTED as the explanation for the observed transfer failure.

## Competing hypotheses

### H1 — Representation / estimation failure

Question:
Does NABU fail because it cannot recover the low-order effects that are actually present?

Tests on low-order fit data:
- exact singleton effect vs learned B2 main effect;
- exact clean-pair e2 vs learned B3 pair contribution.

Report:
- Spearman
- Pearson
- exact std
- learned std
- learned/exact std ratio
- shared identity count

Decision:
- directionally recovered if Spearman >= 0.50 with at least 30 shared identities.
- H1 is SUPPORTED if either singleton or pair direction fails this criterion on both datasets.
- H1 is REJECTED as the primary explanation if both singleton and pair directions pass on both datasets.
- otherwise UNRESOLVED.

Scale shrinkage is reported separately and is not itself classified as representation failure.

### H2 — Coverage / support bias

Question:
Are the apparent higher-order failures mostly caused by which variants have complete low-order evidence?

Within each mutation order with at least 100 supported and 100 unsupported rows:
compare B2 Spearman and target distribution between exact-O2-supported and unsupported rows.

Material support effect:
absolute supported-vs-unsupported B2 Spearman difference >= 0.10.

Decision:
- SUPPORTED if material in the same direction on both datasets.
- REJECTED if <0.10 on both datasets where evaluable.
- otherwise UNRESOLVED.

No across-order mixture comparison may be used to decide H2.

### H3 — Context-free pair portability

Question:
Is an exact pair effect measured in a double reusable as the same additive term in higher-order variants?

Arms:
- exact O1
- exact O2 using measured pair e2
- shuffled-pair O2, seed 161, preserving pair-effect distribution

Evaluate by mutation order on rows with complete exact O2 support.

Evidence against portability:
- exact O2 fails to improve Spearman over O1 by >= 0.05, OR
- exact O2 Spearman is <= shuffled O2 Spearman + 0.02.

Decision:
- SUPPORTED (non-portability) only if the criterion holds on both datasets in at least one matched/evaluable mutation order.
- REJECTED only if exact O2 improves over O1 by >=0.05 and beats shuffled by >0.02 on both datasets in matched/evaluable orders.
- otherwise UNRESOLVED.

### H4 — Context-dependent higher-order interaction

Question:
Do pair effects change materially when a third mutation is present?

For fully O2-supported triples:
e3 = y(abc) - O2(abc)

For each pair in the triple:
conditional_pair = e2(pair) + e3(triple)

Report:
- sign-flip rate
- Spearman(base e2, conditional_pair)
- Pearson(base e2, conditional_pair)
- e3 RMS / target SD

Material context dependence:
- sign-flip rate >= 0.20 AND
- e3 RMS / target SD >= 0.50.

Decision:
- SUPPORTED only if material on both datasets with >=100 supported triples each.
- REJECTED only if below both thresholds on both datasets.
- otherwise UNRESOLVED.

### H5 — Observation-scale / composition nonlinearity

Question:
Is the failure substantially reduced when composition is evaluated on a standard monotonic target scale learned from low-order data only?

Prespecified target transforms fit using fit targets only:
1. identity
2. Yeo-Johnson (standardize=False)
3. robust asinh: asinh((y - median_fit) / MAD_fit), with MAD fallback to std only if MAD==0

For each transform:
- compute transformed WT/single/pair exact components from fit data;
- evaluate transformed O1/O2 on revealed higher-order targets transformed with the already-fit transform;
- no transform parameter may use higher-order targets.

Material improvement over identity:
- O2 Spearman gain >= 0.10 AND
- O2 normalized RMSE reduction >= 20%.

Decision:
- SUPPORTED if at least one non-identity transform meets both criteria on both datasets in a matched/evaluable mutation order.
- REJECTED if neither non-identity transform reaches either criterion on either dataset where evaluable.
- otherwise UNRESOLVED.

This is a falsification diagnostic, not permission to tune a transform on NucB.

### H6 — Low-order information insufficiency / non-identifiability

This hypothesis is intentionally hard to prove.

The audit may only:
- REJECT a strong claim of "no usable structure" if revealed higher-order residuals are predictably structured under strict cross-validation using only features derivable from low-order evidence.
- otherwise mark UNRESOLVED.

It must NOT declare fundamental non-identifiability from failure of one model.

## Global root-cause rule

No single hypothesis may be called ROOT CAUSE unless:
1. Stage-0 invariants pass on both datasets;
2. the hypothesis is SUPPORTED on both datasets with adequate coverage;
3. at least one competing hypothesis is REJECTED by a discriminative test;
4. the result does not depend on mutation-order mixing;
5. the result is reproduced deterministically.

Otherwise final status must be:
ROOT_CAUSE_NOT_YET_IDENTIFIED

## Outputs

- FALSIFICATION_MATRIX.json
- TRPB_INVARIANTS.json
- IRED_INVARIANTS.json
- DATASET_DIAGNOSTICS.json
- RUN_MANIFEST.json
- OUTPUT_HASHES.json

All failed and unresolved results are preserved.
