# NABU External Higher-Order Prior Transfer Audit

Date: 2026-09-28
Branch: exp/nabu-root-cause-audit
Datasets: revealed TrpB + revealed GB1
NucB: forbidden
Phase 3: closed
Architecture selection: none

## Question

WT + singles + doubles do not uniquely identify unrestricted e3+ terms.

Can a SIMPLE UNIVERSAL ALGEBRAIC PRIOR, learned from higher-order observations in one revealed protein landscape, transfer zero-shot to another protein landscape using only the target landscape's low-order evidence?

This is a prior-identification audit, not a NABU replacement model.

## Directional transfer

Two independent directions:

1. GB1 triples -> TrpB triples
2. TrpB triples -> GB1 triples

For each direction:
- source may use source triple labels to learn the external prior;
- target triple labels are NEVER used for fitting, transform selection, hyperparameter selection, feature normalization, or early stopping;
- target labels are used only once for evaluation.

## Scale

Use robust_asinh only.

Its parameters are fit independently from each landscape's LOW-ORDER targets (orders 0-2) only.

No target higher-order labels fit scale parameters.

## Exact low-order evidence

For each landscape on robust_asinh scale:
- derive WT;
- exact singleton e1;
- exact clean-pair e2.

Use only triples with complete e1/e2 support.

## Dimensionless feature vector

Let sigma_low be the standard deviation of robust-asinh LOW-ORDER targets for that landscape.

All magnitude features are divided by sigma_low.

For each supported triple:

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

No sequence identity, protein identity, target triple fitness, or higher-order-derived feature is included.

Target prior quantity:
e3 / sigma_low

## Fixed learner

RandomForestRegressor:
- n_estimators = 300
- min_samples_leaf = 50
- max_features = "sqrt"
- random_state = 161
- n_jobs = 1

No hyperparameter search.

## Negative control

For each transfer direction:
- permute SOURCE e3 targets with seed 161;
- train the identical model;
- evaluate on the same target triples.

## Primary metrics

On target triples:
- Spearman(real transferred prior, true target e3)
- Pearson
- RMSE in dimensionless e3 units
- shuffled-control Spearman
- real minus shuffled Spearman

Also report correlation after removing the target landscape's simple affine e3~O2 component for diagnostic purposes only.

## Gate

UNIVERSAL_ALGEBRAIC_PRIOR_SUPPORTED only if BOTH transfer directions satisfy:
- real target Spearman >= 0.30
- real minus shuffled Spearman >= 0.20

UNIVERSAL_ALGEBRAIC_PRIOR_REJECTED_FOR_THIS_FEATURE_CLASS if BOTH directions have:
- real target Spearman < 0.10
- AND real minus shuffled < 0.10

Otherwise:
UNIVERSAL_ALGEBRAIC_PRIOR_UNRESOLVED

## Interpretation

A rejection applies only to this explicit dimensionless low-order aggregate feature class and fixed learner.

It does NOT reject sequence-aware, structural, or pretrained biological priors.

If rejected, the next justified test is a sequence-aware pretrained prior trained only on low-order labels of the target landscape.

If supported, the learned prior still requires further held-out-landscape validation before architecture use.

## Outputs

- PRIOR_TRANSFER_MATRIX.json
- GB1_TO_TRPB.json
- TRPB_TO_GB1.json
- RUN_MANIFEST.json
- OUTPUT_HASHES.json
