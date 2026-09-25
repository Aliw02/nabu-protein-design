# NABU V9.0 — Pair Transfer Development Protocol

Status: DEVELOPMENT / NOT PHASE 3 / NOT BLIND VALIDATION

Base:
- frozen V8.3 core: c8afdcd7698231d95a39ef13a3fe22b6e3f507c5
- Phase-2 first blind reveal remains immutable.
- V9 development branch starts after Phase-2 closure and does not rewrite any Phase-2 claim.

## Problem isolated by Phase 2D

FLIP2 IRED two-to-many showed:
- 4178 official higher-order test variants;
- 3648 have all main mutation identities represented in fit/train;
- only 131 satisfy frozen V8.3 exact-pair scoreability;
- 4047 therefore abstained under the authoritative Phase-2D execution.

The first V9 experiment changes one mechanism only:
generalization of pair residual evidence to unseen mutation pairs.

AAV indel/insertion semantics are explicitly out of scope for V9.0.
They are a separate future track so the first V9 experiment does not mix two architectural changes.

## V9.0 hypothesis

Frozen V8.3 stores exact pair residual memories.
V9.0 keeps exact pair memories where available, but learns a factorized interaction propensity
for individual mutation identities from cross-fitted double-mutant residuals.

For a previously unseen pair (i,j), the transferred pair residual is estimated as:

    pair_transfer(i,j) = (bias + u_i + u_j) * sqrt(conf_i * conf_j)

where u_i/u_j are minimum-norm least-squares interaction propensities and
conf_i = support_i / (support_i + 1).

Observed exact pairs continue to use the unchanged V8.3 pair memory.
Only missing pair edges use the transferred estimate.

No biochemical word list, dataset-name rule, target-position rule, hand-selected residue class,
or post-hoc threshold is introduced.

## Fit semantics

- Fit data: official IRED train rows with validation != True (3746).
- Held-out validation remains reported separately.
- IRED test labels are already revealed by Phase 2D and may now be used only as a
  post-blind development diagnostic; they no longer count as external validation.
- Frozen V8.3 B2 main memory and exact B3 pair memory are unchanged.
- Factorized pair targets come only from double-mutant fit rows.
- Pair target = label - cross-fitted V8.3 B2 prediction.
- Factorization uses deterministic numpy least squares with no tuned rank or regularization hyperparameter.

## Arms

1. PHASE2_V8_3_ABSTENTION — exact frozen Phase-2 scoreability behavior.
2. B2_MAIN_ONLY — additive frozen V8.3 main effects where all main identities are known.
3. V9_PAIR_TRANSFER — exact V8.3 pairs plus factorized estimates for missing pair edges.
4. V9_PAIR_TRANSFER_PERMUTED_CONTROL — same algorithm with pair residual targets deterministically permuted.

Rows with unseen main mutation identities use the fit global mean in arms 2–4.
No target-aware filtering is allowed.

## Development metrics

Full 4178-row IRED test:
- Spearman
- NDCG
- Top-1% hits / recall / enrichment
- normalized Top-1% regret

Also report:
- all-main coverage;
- exact-pair contribution count;
- transferred-pair contribution count;
- factor-node coverage;
- score distribution diagnostics;
- per-mutation-count metrics.

## Decision rule for continuing V9

V9.0 is considered mechanistically promising only if:
- V9_PAIR_TRANSFER improves full-test Spearman over B2_MAIN_ONLY;
- it improves at least one elite-selection metric without worsening normalized regret;
- it beats the permuted-control arm on Spearman;
- predictions remain finite and deterministic;
- the gain is not confined to the original 131 strict-scoreable rows.

This is a development decision only. It cannot retroactively change Phase 2D.

## Phase boundary

Phase 3 remains unopened.
If V9 development becomes stable, a new untouched external benchmark must be reserved
before any renewed external scientific claim.
