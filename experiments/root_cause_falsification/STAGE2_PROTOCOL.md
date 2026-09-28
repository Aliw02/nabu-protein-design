# NABU Root-Cause Falsification — Stage 2 (GB1 + TrpB)

Date: 2026-09-28
Branch: exp/nabu-root-cause-audit

## Why GB1

IRED has zero higher-order test rows with complete exact-O2 support under its official two-to-many split, so it cannot adjudicate H2-H5 without changing the data contract.

GB1 is already revealed development evidence in this project and has a dense four-site combinatorial landscape. It is used only to test whether the TrpB mechanisms reproduce on a second revealed landscape.

Source:
- repo: rcronk/mutate
- commit: 0861a3f1f6c58a9597fde689eb7f9a07b38e2638
- path: ridge/data/gb1_wu2016_landscape.csv
- Git blob SHA: a1e9f5146a9441210eaf88accfa660f02d80cbb7
- file size: 2,993,024 bytes
- source-documented WT: VDGV
- target: fitness

NucB remains forbidden.
Phase 3 remains closed.

## Diagnostic split

GB1 is not treated as blind.

Mutation order is Hamming distance from WT VDGV.

- fit/low-order evidence: order 0, 1, 2
- test/higher-order evidence: order 3, 4

No higher-order target is used to define exact e1/e2 components or target transforms.

## Frozen thresholds

Stage-1 thresholds are reused unchanged.

### H2 Coverage/support bias
Within each order with >=100 exact-O2-supported and >=100 unsupported rows:
material if |B2 Spearman supported - unsupported| >= 0.10.

Cross-dataset support requires same-direction material effect on TrpB and GB1.

### H3 Context-free pair non-portability
On exact-O2-supported rows:
non-portability evidence if either:
- O2 Spearman - O1 Spearman < 0.05
- O2 Spearman <= shuffled-O2 Spearman + 0.02

Supported only if this holds on TrpB and GB1 in at least one matched order.

### H4 Context-dependent higher-order interaction
On exact-O2-supported triples:
- sign-flip rate >= 0.20
- e3 RMS / target SD >= 0.50
- at least 100 triples

Supported only if material on TrpB and GB1.

### H5 Observation-scale/composition nonlinearity
Reuse exactly:
- identity
- Yeo-Johnson fit on low-order targets only
- robust asinh fit on low-order targets only

Material improvement:
- O2 Spearman gain >= 0.10
- normalized RMSE reduction >= 20%

Supported only if the same prespecified non-identity transform is material on TrpB and GB1 in a matched order.

## H6 Residual-structure diagnostic

This does NOT prove fundamental identifiability.

Question:
After exact O2 is built from low-order evidence, is the higher-order residual predictably structured using only aggregate features derived from low-order evidence?

Rows:
- exact-O2-supported triples and quads only
- evaluated separately by mutation order and pooled

Fixed feature vector (no amino-acid one-hot, no sequence model, no higher-order target-derived feature):
1. mutation order
2. O1
3. O2
4. sum e1
5. mean |e1|
6. max |e1|
7. std e1
8. sum e2
9. mean |e2|
10. max |e2|
11. std e2
12. min e2
13. max e2
14. fraction positive e2

Target:
residual = true higher-order fitness - exact O2

Model:
RandomForestRegressor(
  n_estimators=200,
  min_samples_leaf=20,
  max_features="sqrt",
  random_state=161,
  n_jobs=1
)

CV:
- deterministic 5-fold assignment from FNV-1a(candidate_id + "|NABU_ROOT_CAUSE_H6|161")
- each row predicted only by a model that did not train on that row
- same folds/features/model for shuffled residual control
- shuffled residual uses seed 161

Evidence that residual has usable low-order-derived structure:
- real OOF Spearman >= 0.30
- real minus shuffled OOF Spearman >= 0.20
- both conditions must hold on BOTH TrpB and GB1 pooled supported higher-order rows

If satisfied:
H6_LOW_ORDER_INFORMATION_INSUFFICIENCY = REJECTED AS A STRONG "NO USABLE STRUCTURE" CLAIM.

If not:
H6 remains UNRESOLVED.
It is never marked SUPPORTED by this test.

## Root-cause language rule

Even if H3/H4/H5 reproduce:
- do not call a single ROOT CAUSE unless the cross-dataset discrimination rules are satisfied.
- if multiple mechanisms remain compatible, output ROOT_CAUSE_NOT_YET_IDENTIFIED.

## Outputs

- STAGE2_MATRIX.json
- GB1_INVARIANTS.json
- GB1_DIAGNOSTICS.json
- H6_RESIDUAL_STRUCTURE.json
- RUN_MANIFEST.json
- OUTPUT_HASHES.json
