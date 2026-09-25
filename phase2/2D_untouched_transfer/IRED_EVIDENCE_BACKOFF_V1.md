# Phase 2D-B IRED Full-Test Handling Freeze

**Version:** `IRED_EVIDENCE_BACKOFF_V1`  
**Frozen before any IRED target reveal.**

## Trigger

The label-blind structural preflight on the official 4,178-row IRED test set found:

- frozen V8.3 scoreable: 131 / 4,178;
- all main mutation identities supported: 3,648 / 4,178;
- at least one supported internal pair: 133 / 4,178.

The preregistered protocol forbids scoreable-only substitution for the full-test
primary metrics, so a complete-test handling rule is required before labels are
read.

## Frozen evidence backoff

For every official test identity:

1. **ROUTER**
   - if all mutation main effects are supported and at least one internal pair
     is supported;
   - use the unchanged frozen `V8_3_ADAPTIVE_ROUTER` score.

2. **B2_ADDITIVE**
   - if all mutation main effects are supported but no internal pair is
     supported;
   - use the unchanged frozen B2 additive score.

3. **GLOBAL_MEAN_ABSTENTION**
   - if one or more mutation main effects are unseen in fit/train;
   - use the fit/train global mean as a neutral abstention score.

No coefficients, thresholds, memories, support requirements, or router logic are
changed. The policy selects the deepest already-frozen evidence level whose
structural prerequisites are satisfied.

## Full-test metrics

Primary Spearman and NDCG and all frozen secondary metrics are computed on all
4,178 official test identities using the resulting finite prediction vector.

The original V8.3 strict-scoreability subset is reported only as a coverage and
diagnostic result and cannot replace the full-test primary metrics.

## Tie breaking

Stable source order is used for equal prediction values. No target-aware
tie-breaking is permitted.

## Scientific interpretation

This rule is an evaluation-time abstention/backoff policy required to make the
already-frozen hierarchy total over the official IRED test domain. It is not a
new learned model and is frozen before any IRED target value is read.
