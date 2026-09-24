# NABU V8 Adaptive Higher-Order Development Protocol

Status: DEVELOPMENT ONLY — NOT A CANONICAL FREEZE.

## Purpose

V7 showed that most of the gain over additive scoring came from the raw pairwise residual model, while the consensus layer only made a small ranking change. V8 therefore tests a different mechanism: whether residual signal remains after pairwise modeling and can be recovered from supported higher-order mutation subsets.

PHOT_CHLRE_Chen_2023 is now a development/debug assay because its outcomes have already been revealed. No result from this branch should be described as a fresh blind confirmation.

## Hierarchy

1. B2_ADDITIVE
   - Global visible mean.
   - Mutation main effects with n/(n+1) shrinkage.

2. B3_RAW_PAIR
   - Pair residual memory fit on y - B2.
   - Same raw-pair summation semantics used by V7 for a direct comparison.

3. B4_TRIPLET_RESIDUAL
   - Fit triplet memory on residual y - B3.
   - Require support >= 2.
   - Shrinkage confidence n/(n+2).
   - Candidate contribution is an adaptive confidence-weighted residual averaged across supported triplets and scaled by supported-subset coverage.

4. B5_ADAPTIVE_HIGHER_ORDER
   - Fit quartet memory on residual y - B4.
   - Require support >= 2.
   - Shrinkage confidence n/(n+4).
   - Candidate contribution uses the same adaptive support/coverage rule.
   - No exact full-sequence memorization term is used. With 3-5 mutation candidates, the highest generalizable proper-subset order is quartet.

## Controls

Two true permutation controls replace the V7 row-rotation control:

- GLOBAL_PERMUTED_B5: deterministic full permutation of visible labels with seed 161.
- COUNT_STRATIFIED_PERMUTED_B5: deterministic permutation within mutation-count strata, preserving the 3/4/5-mutation count distribution while destroying mutation identity-to-fitness mapping.

## Outputs

- V8_DEV_RESULTS.json
- V8_DEV_REVEALED_PREDICTIONS.csv
- V8_DEV_TOP50.csv
- V8_DEV_MEMORY_DIAGNOSTICS.json

## Breakpoint rule

Do not freeze this branch merely because a metric improves. It becomes a breakpoint candidate only if the higher-order layer gives a useful and reproducible gain over B3_RAW_PAIR without relying on broken controls or trivial candidate overlap. If that happens, the architecture is reimplemented as a sealed two-stage run and tested on a fresh unrevealed landscape.
