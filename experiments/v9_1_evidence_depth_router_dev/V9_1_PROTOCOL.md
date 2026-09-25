# NABU V9.1 — Graceful Evidence-Depth Router

Status: DEVELOPMENT / NOT PHASE 3 / NOT BLIND VALIDATION

Base:
- frozen V8.3 core remains unchanged;
- V9.0 factorized unseen-pair transfer is rejected and preserved;
- IRED is development-only because Phase 2D already revealed its targets.

## Problem

Phase 2D showed that frozen V8.3 hard scoreability collapses on IRED:
131 / 4178 test variants are strict-scoreable even though 3648 / 4178 have all
main mutation identities represented.

V9.0 showed that inventing unseen pair residuals through node factorization
degrades the stronger B2 main-effect signal.

## Single V9.1 change

Replace the hard exact-pair scoreability requirement with an evidence-depth
router.

For a candidate:

1. If any main mutation identity is unseen:
   use the fit global mean as neutral abstention.
2. If all main identities are known:
   B2 additive is always defined.
3. Exact B3 pair correction is allowed only if fit-only out-of-fold evidence
   shows B3 ranks visible fit rows better than B2:
      Spearman(OOF_B3, fit_target) > Spearman(OOF_B2, fit_target)
4. If that OOF gate is closed, all all-main candidates remain at B2.
5. No unseen pair effect is estimated.
6. No dataset-name rule, biochemical word list, position rule, or test-target
   threshold is introduced.

For all-main candidates with no exact pair memory, B3 numerically reduces to B2.

## Arms

- PHASE2_V8_3_ABSTENTION
- B2_MAIN_ONLY
- B3_AVAILABLE_ALWAYS (mechanism control)
- V9_1_EVIDENCE_DEPTH_ROUTER

## Development metrics

Full 4178-row IRED diagnostic:
- Spearman
- NDCG
- Top-1% hits / recall / enrichment
- normalized Top-1% regret

Diagnostics:
- OOF B2 and B3 Spearman on fit;
- pair-depth gate state;
- all-main coverage;
- strict frozen scoreability;
- strict-subset B2/B3 comparison;
- non-strict all-main comparison;
- held-out validation metrics;
- prediction hashes and deterministic replay.

## Development decision

V9.1 is retained as the next architecture base only if:

- predictions are finite and deterministic;
- all-main coverage expands scoring beyond frozen strict coverage;
- full-test Spearman is >= B2_MAIN_ONLY;
- full-test Top-1% hits are >= B2_MAIN_ONLY;
- normalized Top-1% regret is <= B2_MAIN_ONLY.

If the OOF gate opens, the result must also report whether the exact-pair
correction helps or hurts the strict subset; that diagnostic does not retune
the gate after test observation.

Equality with B2 is acceptable for retention because the architectural goal is
graceful evidence coverage without inventing unsupported interactions.
A fresh untouched dataset is required before any external scientific claim.

## Phase boundary

Phase 3 remains unopened.
