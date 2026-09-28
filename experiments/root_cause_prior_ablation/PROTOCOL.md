# NABU Prior Transfer Final Ablation

Date: 2026-09-28
Branch: exp/nabu-root-cause-audit
NucB: forbidden
Phase 3: closed

## Purpose

The cross-landscape triple prior transferred strongly:
- GB1 -> TrpB
- TrpB -> GB1

Before any architecture is allowed, test whether that result is merely a scalar/nonlinear O2 cancellation shortcut or whether low-order interaction-shape features carry transferable information beyond O2.

This is the stopping test for the current root-cause investigation.

## Frozen data/scale/model

Use the exact same:
- TrpB and GB1 sources;
- robust_asinh fit from target/source low-order labels only;
- exact-O2-supported triples;
- dimensionless normalization;
- RandomForestRegressor parameters;
- seed 161;
- cross-landscape directions;
as the External Higher-Order Prior Transfer audit.

No hyperparameter changes.

## Feature arms

### O2_ONLY
One feature:
- normalized O2

### SHAPE_ONLY
Exclude O1, O2, sum(e1), sum(e2).

Use:
- mean(abs(e1))
- max(abs(e1))
- std(e1)
- min(e1)
- max(e1)
- mean(abs(e2))
- max(abs(e2))
- std(e2)
- min(e2)
- max(e2)
- fraction(e2 > 0)

All magnitudes remain normalized by source/target low-order sigma as already frozen.

### FULL
The exact 15-feature vector from the prior-transfer audit.

## Evaluation

For each direction and arm:
- target Spearman
- target Pearson
- RMSE
- conditional residual Spearman after independently removing the target O2 affine component from true e3 and prediction.

No target labels are used for model fitting.

The residualization is evaluation-only.

## Primary stopping gate

TRANSFERABLE_INTERACTION_SHAPE_SUPPORTED only if BOTH directions satisfy:

1. FULL conditional residual Spearman >= 0.30
2. FULL conditional residual Spearman - O2_ONLY conditional residual Spearman >= 0.10
3. SHAPE_ONLY target Spearman >= 0.20

If this passes:
- stop adding root-cause diagnostics;
- the evidence supports a transferable low-order interaction-shape prior beyond scalar O2 cancellation;
- the next work is one architecture derived from that prior, tested on revealed development only before NucB.

If it fails:
- do not build the algebraic prior architecture;
- classify the previous transfer as predominantly cancellation/scale shortcut;
- the next prior class to test is sequence-aware pretrained context.

## Outputs

- PRIOR_ABLATION_MATRIX.json
- GB1_TO_TRPB_ABLATION.json
- TRPB_TO_GB1_ABLATION.json
- RUN_MANIFEST.json
- OUTPUT_HASHES.json
