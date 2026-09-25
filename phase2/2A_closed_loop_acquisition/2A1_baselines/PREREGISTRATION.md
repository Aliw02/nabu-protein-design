# Phase 2A.1 Baseline Tournament — Preregistration

## Purpose

Characterize fair acquisition controls before designing a new Phase-2 controller.

This stage is not allowed to tune a new controller.

## Landscape

Initial development landscape:

- RhlA training pool only;
- full pool used as a virtual assay oracle;
- no RhlA hidden-triple test truth is used for acquisition.

## Shared campaign protocol

All baselines use:

- identical identity-only 5% bootstrap;
- identical frozen V8.3 core;
- identical candidate pool;
- identical oracle;
- identical batch size;
- identical target measurement count;
- deterministic tie-breaking;
- identical evaluation checkpoints.

Default tournament:

- seed fraction: 5%;
- target fraction: 20%;
- batch size: 16.

## Frozen baselines

1. **deterministic_random**
   - fixed pseudorandom candidate permutation;
   - independent of labels and model predictions.

2. **greedy_v83**
   - refits after every batch;
   - selects by the current within-pool V8.3 ordering.

3. **exploration_only**
   - uses the historical structural support deficit;
   - main and pair support only;
   - no exploitation score.

4. **historical_50_50**
   - historical structural support-deficit percentile rank;
   - current V8.3 exploitation percentile rank;
   - fixed 0.5 / 0.5 rank fusion.

5. **static_initial_v83**
   - freezes the initial post-bootstrap V8.3 ranking;
   - later refits occur for diagnostics but do not change acquisition order.

## Score-scale rule

V8.3 router output is treated as an ordering, not a globally calibrated fitness value.

Policies use within-pool ordering/rank semantics.

## Primary metrics

- best true fitness discovered vs measurements spent;
- best-so-far regret;
- cumulative Top-1% hits;
- Top-1% hits acquired after bootstrap;
- first Top-1% measurement checkpoint;
- normalized post-bootstrap discovery AUC;
- diversity among discovered Top-1% variants.

Diagnostics:

- scoreable coverage;
- visible B3 OOF Spearman;
- router trajectory.

## Interpretation boundary

This single RhlA tournament establishes controls.

It does not select the final Phase-2 controller and does not support a universal active-learning claim.

No new controller weight or rule may be tuned from the final RhlA ordering and then reported as untouched evidence.
