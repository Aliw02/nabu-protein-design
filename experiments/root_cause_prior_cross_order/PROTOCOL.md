# NABU Cross-Landscape / Cross-Order Prior Transfer Audit

Date: 2026-09-28
Branch: exp/nabu-root-cause-audit
Datasets: revealed TrpB + revealed GB1
NucB: FORBIDDEN
Phase 3: CLOSED
Architecture selection: NONE

## Question

The previous audit showed that a fixed dimensionless algebraic prior trained on source-protein triples transfers strongly to target-protein triples, including after removing a simple e3~O2 affine component.

This audit asks the harder discriminative question:

> Can a prior learned ONLY from triples in protein A predict the missing higher-order residual of QUADS in protein B, using only protein B's WT/single/double evidence?

No target triple or quad labels are used for training, feature normalization, transform fitting, hyperparameter selection, or early stopping.

## Directions

1. GB1 triples -> TrpB quads
2. TrpB triples -> GB1 quads

Triple-to-triple transfer is recomputed only as a replication check and cannot determine the primary decision.

## Scale

robust_asinh only.

Each landscape fits center/scale from its own LOW-ORDER targets (orders 0-2) only.

## Features

Exactly the prior-transfer feature class, dimensionless by sigma_low:

1. O1 / sigma_low
2. O2 / sigma_low
3. sum(e1) / sigma_low
4. mean(abs(e1)) / sigma_low
5. max(abs(e1)) / sigma_low
6. std(e1) / sigma_low
7. min(e1) / sigma_low
8. max(e1) / sigma_low
9. sum(e2) / sigma_low
10. mean(abs(e2)) / sigma_low
11. max(abs(e2)) / sigma_low
12. std(e2) / sigma_low
13. min(e2) / sigma_low
14. max(e2) / sigma_low
15. fraction(e2 > 0)

No sequence identity, protein identity, target higher-order label, or target higher-order-derived feature.

Target residual:
higher_order_residual / sigma_low = (transformed target - exact O2) / sigma_low

## Fixed learner

RandomForestRegressor:
- n_estimators=300
- min_samples_leaf=50
- max_features="sqrt"
- random_state=161
- n_jobs=1

Train ONLY on source-protein supported triples.

## Controls

### Shuffled-source residual control
Permute source triple residual labels with seed 161 and train identical RF.

### Source affine cancellation control
Fit on source triples only:

residual = a + b * (O2 / sigma_low)

Freeze a,b and apply to target quads.

This control directly tests whether any apparent cross-order transfer is only the simple cancellation law already observed.

### Target B2 reference
Fit frozen NABU B2 from target WT/single/double evidence only and score target quads.

This is a reference, not training for the transferred prior.

## Metrics on target quads

Residual prediction:
- Spearman(real prior, true target residual)
- Pearson
- RMSE
- shuffled-control Spearman
- source-affine-control Spearman
- real minus shuffled
- real minus source-affine

Reconstructed fitness:
predicted transformed fitness = exact O2 + predicted residual

Report:
- Spearman vs true target fitness
- exact O2 Spearman
- target B2 Spearman
- shuffled-prior reconstructed Spearman
- source-affine reconstructed Spearman

## Predeclared primary gate

CROSS_LANDSCAPE_CROSS_ORDER_PRIOR_SUPPORTED only if BOTH transfer directions satisfy ALL:

1. quad residual Spearman >= 0.30
2. real minus shuffled residual Spearman >= 0.20
3. real minus source-affine residual Spearman >= 0.10
4. reconstructed quad fitness Spearman >= exact O2 Spearman + 0.10
5. reconstructed quad fitness Spearman >= target B2 Spearman - 0.02

CROSS_LANDSCAPE_CROSS_ORDER_PRIOR_REJECTED_FOR_THIS_FEATURE_CLASS only if BOTH directions satisfy:
- residual Spearman < 0.10
- AND real minus shuffled < 0.10

Otherwise:
CROSS_LANDSCAPE_CROSS_ORDER_PRIOR_UNRESOLVED

## Interpretation

If supported:
- target low-order evidence alone still does not identify unrestricted e3/e4;
- but a transferable external higher-order prior can supply useful missing structure across BOTH protein identity and mutation order;
- this justifies testing a frozen prior-augmented architecture on another revealed development landscape before NucB.

If unresolved/rejected:
- do NOT patch weights;
- move to a sequence-aware/pretrained prior or additional higher-order measurements.

This audit cannot open Phase 3 and cannot consume NucB.

## Outputs

- CROSS_ORDER_MATRIX.json
- GB1_TRIPLES_TO_TRPB_QUADS.json
- TRPB_TRIPLES_TO_GB1_QUADS.json
- RUN_MANIFEST.json
- OUTPUT_HASHES.json
