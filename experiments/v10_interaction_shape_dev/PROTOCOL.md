# NABU V10.0 — Interaction-Shape Prior Development Protocol

Date: 2026-09-28
Branch: exp/nabu-v10-interaction-shape-dev
Status: development only
NucB: FORBIDDEN
Phase 3: CLOSED

## Purpose

Build exactly ONE architecture from the preregistered final-ablation result:

TRANSFERABLE_INTERACTION_SHAPE_SUPPORTED

This is not another root-cause diagnostic and not a V9 patch.

If V10.0 fails the frozen development gate, V10.0 is rejected. No V10.1 weight tuning, threshold tuning, blend tuning, order-specific patch, or local reranker is allowed from the same revealed results.

## Architectural principle

B2 remains the global ranking backbone.

The external prior is allowed to predict only the RANK ERROR of B2 using the interaction-shape features that survived the final ablation.

For a source landscape:

1. fit B2 from source WT + singles + doubles only;
2. score every higher-order source candidate with B2;
3. convert source B2 scores to global percentile ranks;
4. convert source true higher-order fitness to global percentile ranks;
5. define source rank correction target:

   rank_residual = true_percentile_rank - B2_percentile_rank

6. train a frozen interaction-shape prior on supported source higher-order candidates.

For a target landscape:

1. fit B2 from target WT + singles + doubles only;
2. compute target B2 percentile rank over all higher-order target candidates;
3. compute interaction-shape features from target LOW-ORDER evidence only;
4. predict rank correction on exact-pair-supported target candidates;
5. unsupported candidates receive correction exactly 0;
6. final V10 score:

   V10_score = B2_percentile_rank + predicted_rank_correction

No clipping, learned blending coefficient, hand-tuned weight, or target-label calibration.

## Cross-landscape development directions

Two independent directions:

1. GB1 revealed higher-order -> TrpB target
2. TrpB revealed higher-order -> GB1 target

In each direction:
- source higher-order labels may train the prior;
- target higher-order labels are NEVER used for fitting, feature normalization, transform fitting, hyperparameter selection, early stopping, or correction calibration;
- target labels are used once for evaluation.

## Low-order scale

Use robust_asinh exactly as frozen in the final ablation.

Parameters are fit independently from each landscape's low-order targets only:

z = asinh((y - median_low_order) / MAD_low_order)

MAD fallback to low-order std only if MAD == 0.

Let sigma_low be std(z_low_order).

## V10 interaction-shape features

Use EXACTLY the SHAPE_ONLY feature set that passed the final ablation.

All magnitude features are divided by sigma_low.

1. mean(abs(e1))
2. max(abs(e1))
3. std(e1)
4. min(e1)
5. max(e1)
6. mean(abs(e2))
7. max(abs(e2))
8. std(e2)
9. min(e2)
10. max(e2)
11. fraction(e2 > 0)

No O1.
No O2.
No sum(e1).
No sum(e2).
No sequence identity.
No protein identity.
No target higher-order-derived feature.
No mutation-order-specific model.

A candidate is V10-supported only when every required singleton and pair term exists in target low-order evidence.

## Fixed prior learner

RandomForestRegressor:
- n_estimators = 300
- min_samples_leaf = 50
- max_features = "sqrt"
- random_state = 161
- n_jobs = 1

No hyperparameter search.

## Controls

### B2
Unmodified B2 backbone.

### O2_ONLY correction control
Same rank-residual target and same RF, but only normalized O2 is supplied as the prior feature.

This tests whether V10 is merely reintroducing the cancellation shortcut rejected by the final ablation.

### Shuffled-rank-residual control
Permute source rank-residual labels with seed 161 and train the V10 feature model unchanged.

## Metrics

Report on the complete target higher-order set:
- Spearman
- NDCG
- Top-1% hits
- normalized regret
- best true selected

Report separately for mutation order 3 and 4:
- Spearman
- NDCG
- Top-1% hits

Also report:
- exact V10 support count/fraction;
- supported-region B2 vs V10;
- unsupported-region B2 vs V10 (must be identical);
- O2_ONLY correction control;
- shuffled-rank-residual control.

## Determinism

Train V10 twice from identical inputs in each direction.

Maximum absolute difference between replayed V10 prediction vectors must be 0.0.

## Frozen development gate

V10_DEV_PASS only if BOTH cross-landscape directions satisfy ALL:

1. whole-target V10 Spearman > B2 Spearman;
2. whole-target V10 NDCG >= B2 NDCG;
3. whole-target V10 Top-1% hits >= B2 Top-1% hits;
4. order-3 V10 Spearman >= order-3 B2 Spearman;
5. order-4 V10 Spearman >= order-4 B2 Spearman;
6. whole-target V10 Spearman > O2_ONLY correction-control Spearman;
7. whole-target V10 Spearman > shuffled-control Spearman;
8. unsupported-region V10 scores are exactly equal to B2 percentile scores;
9. deterministic replay max absolute difference == 0.0.

No tolerance may be added after observing results.

## Decision

If every condition passes in BOTH directions:

V10_DEV_PASS

Then:
- freeze V10.0 architecture;
- train the final external prior on BOTH revealed development landscapes;
- save its hashes and protocol;
- DO NOT consume NucB in this run;
- NucB blind becomes the next separate scientific action.

If any condition fails:

V10_DEV_REJECT

Then:
- reject V10.0;
- do not patch it on revealed target results;
- do not create V10.1 from the same failure;
- NucB remains untouched.

## Outputs

- V10_DEV_MATRIX.json
- GB1_TO_TRPB_V10.json
- TRPB_TO_GB1_V10.json
- FINAL_PRIOR_FREEZE.json only if V10_DEV_PASS
- FINAL_PRIOR_MODEL.joblib only if V10_DEV_PASS
- RUN_MANIFEST.json
- OUTPUT_HASHES.json
