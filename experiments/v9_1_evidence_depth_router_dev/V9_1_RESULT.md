# NABU V9.1 Evidence-Depth Router — Development Result

Status: RETAIN AS NEXT ARCHITECTURE BASE

Phase 3 remains unopened. Frozen V8.3 remains unchanged.

## Core finding

The fit-only OOF depth gate observed:

- OOF B2 Spearman: 0.32860921025867007
- OOF B3 Spearman: 0.32860921025867007
- delta: 0.0
- pair depth enabled: false

Therefore V9.1 correctly refused to add exact-pair correction and used B2 for
all candidates whose main mutation identities are known.

This is especially important on IRED because every frozen fit pair memory has
support 1. Under cross-fitting, an exact held-out pair is not reusable, so B3
does not provide OOF evidence beyond B2.

## Coverage

- frozen strict V8.3 scoreable: 131 / 4178
- all-main-supported: 3648 / 4178
- all-main without exact pair: 3517
- unseen-main abstentions remaining: 530

V9.1 therefore expands evidence-backed scoring coverage from 3.135% strict
coverage to 87.315% all-main coverage without inventing unseen pair effects.

## Full IRED development diagnostic

| Arm | Spearman | NDCG | Top-1% hits | Normalized regret |
| --- | ---: | ---: | ---: | ---: |
| Phase-2 V8.3 abstention | -0.021771 | 0.942698 | 0 / 42 | 0.231522 |
| B2 main only | 0.169971 | 0.956090 | 3 / 42 | 0.064702 |
| B3 available always | 0.170698 | 0.955141 | 3 / 42 | 0.064702 |
| V9.1 evidence-depth router | 0.169971 | 0.956090 | 3 / 42 | 0.064702 |

V9.1 exactly matches the stronger B2 baseline because the fit-only OOF gate
does not justify B3.

## Strict exact-pair diagnostic

131 rows:

- B2 Spearman: 0.111039
- B3 Spearman: 0.119388

B3 is slightly better after reveal on this tiny subset, but V9.1 does not use
that post-reveal fact to override the fit-only OOF gate.

## Held-out validation diagnostic

662 rows:

- all-main-supported: 613
- exact-pair-supported: 0
- B2 / B3 / V9.1 Spearman: 0.321370

The validation domain independently supports the graceful B2 fallback behavior.

## Decision

V9.1 passes its preregistered development retention rule:

- finite predictions: PASS
- deterministic replay: PASS
- coverage expands beyond strict V8.3: PASS
- Spearman not worse than B2: PASS
- Top-1% hits not worse than B2: PASS
- normalized regret not worse than B2: PASS

Decision:
RETAIN_V9_1_AS_NEXT_ARCHITECTURE_BASE

This is development evidence on already-revealed IRED, not a renewed external
scientific validation. A new untouched dataset is required before any external
claim.

Remaining architectural gap:
530 IRED test variants still contain at least one unseen main mutation identity.
No V9.1 rule guesses those effects.

CI run:
- run: 36196149069
- artifact: 10889927267
- artifact digest: sha256:561411c78357c1f8c29f32d7d0d96e45cd4c69ab7824818735600a93872691e1
