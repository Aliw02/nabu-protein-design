# NABU V9.0 Pair Transfer — Development Result

Status: REJECTED DEVELOPMENT HYPOTHESIS

Phase 3 remains unopened. Frozen V8.3 is unchanged.

## Hypothesis tested

V9.0 attempted to transfer double-mutant residual information to unseen mutation
pairs using a factorized node model:

    pair(i,j) = (bias + u_i + u_j) * sqrt(conf_i * conf_j)

Exact frozen V8.3 pair memories retained priority.

## Data role

FLIP2 IRED two-to-many is no longer blind after Phase 2D.
It is used here only as a development diagnostic.

Split:
- fit: 3746
- validation: 662
- test diagnostic: 4178

## Coverage

- all-main-supported test variants: 3648
- frozen strict V8.3 scoreable: 131
- known-main / non-strict transfer region: 3517
- exact pair contributions: 137
- transferred pair contributions: 20831
- unsupported pair contributions: 5291

## Full-test development metrics

| Arm | Spearman | NDCG | Top-1% hits | Normalized regret |
| --- | ---: | ---: | ---: | ---: |
| Phase-2 V8.3 abstention | -0.021771 | 0.942698 | 0 / 42 | 0.231522 |
| B2 main only | 0.169971 | 0.956090 | 3 / 42 | 0.064702 |
| V9.0 pair transfer | 0.095192 | 0.947637 | 0 / 42 | 0.215702 |
| Permuted pair-transfer control | 0.086144 | 0.944080 | 0 / 42 | 0.217646 |

V9.0 beat its permuted control slightly on Spearman, but it lost clearly to the
simpler B2 main-only arm on every primary discovery signal.

## Non-strict region

On the 3517 all-main-supported rows that failed frozen exact-pair scoreability:

- B2 Spearman: 0.185930
- V9.0 pair-transfer Spearman: 0.101495

The transfer mechanism therefore degraded the exact region it was designed to
repair.

## Validation diagnostic

Held-out validation showed the same direction:

- B2 Spearman: 0.321370
- V9.0 pair-transfer Spearman: 0.260745

## Factor identifiability

- double-mutant rows: 2708
- factor nodes: 1298
- parameters including bias: 1299
- design rank: 1290

The least-squares node factorization is underdetermined. Independent CI runs
preserved identical rankings and scientific metrics but showed tiny raw-float
differences. Numerical stabilization did not make the factor predictions
bit-identical across runners.

## Decision

V9.0 factorized unseen-pair transfer is not promoted.

The strongest development signal is instead that B2 main effects generalize
substantially better than hard abstention, while guessed unseen-pair residuals
destroy part of that signal.

Next architecture direction:
V9.1 graceful evidence-depth routing — use B2 whenever all main identities are
known, and permit deeper exact-pair evidence only through a fit/OOF-only gate.
No unseen pair effect will be invented.

This result does not alter the frozen Phase-2 verdict.
