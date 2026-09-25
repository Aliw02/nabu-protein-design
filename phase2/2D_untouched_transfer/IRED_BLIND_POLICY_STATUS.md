# Phase 2D IRED Blind Policy Status

## Authoritative full-test policy

**IRED_ABSTENTION_V2**

This policy was frozen before target reveal.

- Frozen V8.3 scoreable test rows: use unchanged `V8_3_ADAPTIVE_ROUTER`.
- Frozen V8.3 unscoreable test rows: retain abstention semantics and assign one
  tied ranking score below all finite scoreable router scores.
- No B2 fallback is used.
- No test rows are dropped.
- Primary Spearman and NDCG use all 4,178 official test rows.

## Identity-only alternative not selected

An identity-only implementation of a B2/global-mean evidence backoff was committed
during protocol development. It was rejected before target reveal because it
would change the frozen V8.3 scoreability semantics.

Its code/history remains preserved for audit, but it is not authoritative and
must not be used to report the Phase-2D blind result.

Authoritative runner:

`run_ired_blind_v2.py`
