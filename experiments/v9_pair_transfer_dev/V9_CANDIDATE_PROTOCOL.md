# NABU V9 Candidate Architecture — Clean Calibrated Pair Expansion

Status: DEVELOPMENT CANDIDATE / POST-PHASE-2
Phase 3: CLOSED
Scientific claim: false

## Purpose

Run one decisive higher-order development comparison after the clean-pair and
nested-OOF calibration diagnostics succeeded.

This candidate does not alter frozen V8.3 or the Phase-2 verdict.

## Frozen evidence decomposition

For the official IRED fit/train partition only:

1. WT target:
       y_WT

2. Clean main effect for mutation m when its singleton is measured:
       main(m) = y(single_m) - y_WT

3. Clean observed pair epistasis when both singleton counterparts and the
   double mutant are measured:
       epsilon(i,j) = y(i,j) - y(single_i) - y(single_j) + y_WT

4. Unseen-pair relation:
       raw_epsilon_hat(i,j) = bias + u_i + u_j

5. Calibration:
       calibrated_epsilon_hat = alpha + beta * raw_epsilon_hat

The global alpha/beta are learned only from deterministic outer-fold OOF
predictions over the 1,759 clean fit/train pair rows. No IRED validation or
test target is used to fit the relation or calibration.

## V9 higher-order score

For a candidate with mutations M:

If every mutation has a measured singleton:

    score =
        y_WT
        + sum(main(m) for m in M)
        + sum(pair_effect(i,j) for all unordered pairs in M)

where pair_effect is:

- exact clean measured epsilon(i,j), if that clean pair was observed in fit/train;
- otherwise calibrated predicted epsilon(i,j), if both mutation identities are
  represented in the fitted clean-pair relation;
- otherwise 0.0.

This is the standard second-order expansion. No degree normalization, clipping,
pair-count scaling, learned threshold, mutation-count rule, or dataset-specific
scalar is added.

If any mutation lacks a measured singleton, the entire candidate falls back to
the unchanged V8.3 B2 main-only score. This fallback is identity/evidence based,
is reported explicitly, and prevents silent subsetting.

## Development comparison

Evaluate on the already-revealed IRED test partition as development data only.

Arms:

1. PHASE2_V8_3_ABSTENTION
2. B2_MAIN_ONLY
3. V9_0_RAW_PAIR_TRANSFER
4. V9_1_DEGREE_NORMALIZED_TRANSFER
5. V9_CLEAN_CALIBRATED_PAIR_EXPANSION

Report:
- full-test Spearman;
- NDCG;
- Top-1% hits / recall / enrichment;
- normalized Top-1% regret;
- metrics by mutation count;
- clean-main coverage;
- exact clean-pair contributions;
- calibrated unseen-pair contributions;
- unresolved pair contributions;
- B2 fallback count;
- prediction scale and deterministic replay hash.

## Decisive development rule

The V9 candidate is considered ready to retain as the V9 development architecture
only if ALL are true:

1. full-test Spearman > B2_MAIN_ONLY;
2. full-test Spearman > V9_1_DEGREE_NORMALIZED_TRANSFER;
3. Top-1% hits >= B2_MAIN_ONLY;
4. normalized Top-1% regret <= B2_MAIN_ONLY;
5. on candidates with complete clean-main support, Spearman > B2 on the same rows;
6. Spearman is positive for mutation counts 3, 4, and 5 separately;
7. prediction standard deviation is finite and <= 3x target standard deviation;
8. deterministic replay preserves identical rank order.

If this fails, stop identity-only pair-expansion development and move to a
richer relational sequence/structure representation.

## Boundary

IRED test labels were already revealed in Phase 2D, so this run is development
evidence only. A success cannot repair Phase 2 and cannot open Phase 3.
A new untouched benchmark is required afterward.
