# NABU V9 Contextual ESM Residual — Result

Status: COMPLETE
Strict decision: V9_CONTEXTUAL_ESM_RESIDUAL_REJECT
Phase 3: CLOSED
Scientific claim: false

## Execution

GitHub Actions run:
- 36237040369

Workflow head:
- da557a957209e1c46f3f82188cd1e49302106d08

Artifact:
- 10904493017
- sha256:a52c3338b8eaa01abe463484cb8bd376eaf3a7d2c70b80b0ab05919c62421fcc

All contracts and the workflow completed successfully.

## Context encoder

Frozen encoder:
- facebook/esm2_t6_8M_UR50D
- revision c731040fcd8d73dceaa04b0a8e6329b345b0f5df
- feature dimension 640
- no encoder fine-tuning

Feature:
- global mean final-hidden delta vs WT
- mean mutation-site final-hidden delta vs WT

Context readout:
- RidgeCV
- selected alpha: 1000.0

## Residual-training diagnostic

Training evidence:
- WT: 1
- singles: 1,037
- clean doubles: 1,759
- total residual-training rows: 2,797
- excluded doubles missing singleton counterpart: 949

OOF residual prediction:
- all rows Spearman: 0.1988236587
- all rows RMSE: 1.0077467771
- clean doubles Spearman: 0.2598502855
- clean doubles RMSE: 1.2478930616

The frozen contextual representation therefore contains learnable interaction signal before the higher-order test.

## Full IRED development test

| Arm | Spearman | NDCG | Top-1% hits | Normalized regret |
|---|---:|---:|---:|---:|
| Phase-2 V8.3 abstention | -0.0217706 | 0.9426984 | 0/42 | 0.2315224 |
| B2 main-only | 0.1699714 | 0.9560896 | 3/42 | 0.0647016 |
| V9.1 degree-normalized | 0.1370972 | 0.9528017 | 3/42 | 0.0647016 |
| V9 clean calibrated pair expansion | 0.1551567 | 0.9498469 | 0/42 | 0.2103395 |
| V9 contextual ESM residual | 0.1734523 | 0.9564649 | 1/42 | 0.0647016 |

Important:
- contextual V9 is the first V9 arm to beat B2 on full-test Spearman;
- NDCG is also slightly above B2;
- it retains the same best true selected variant as B2;
- Top-1% overlap is lower: 1 hit versus B2's 3.

## Complete-clean-main region

Rows: 1,871

| Arm | Spearman | NDCG | Top-1% hits | Normalized regret |
|---|---:|---:|---:|---:|
| B2 | 0.2114864 | 0.9480111 | 3/19 | 0.0593089 |
| Clean calibrated pair expansion | 0.1789433 | 0.9338693 | 0/19 | 0.2531222 |
| Contextual ESM residual | 0.2572264 | 0.9504725 | 1/19 | 0.0593089 |

On the rows where the new mechanism actually operates, contextual V9 improves
Spearman over B2 by about 0.04574.

## Mutation-order behavior

Contextual V9 vs B2 Spearman:

- order 3: 0.2282967 vs 0.2300949
- order 4: 0.1747792 vs 0.1779972
- order 5: 0.0327732 vs 0.0267052
- order 6: 0.0710410 vs 0.0443677
- order 7: 0.1504634 vs 0.0626766
- order 8: -0.0854653 vs -0.1628744
- order 9: 0.4181818 vs 0.4666667 (n=10)
- order 10: 0.2581989 vs -0.2581989 (n=4)

The preregistered positive-order-3/4/5 condition passed.

## Elite-selection diagnostic

The strict failure is localized to Top-1% overlap.

Among the 42 true Top-1% variants:
- 19 have complete clean-single support;
- 23 lie in the B2-fallback region.

B2's three Top-1% hits were all complete-clean-main variants with targets:
- 1.4836020520
- 2.2558416679
- 2.4979888214

Contextual V9 retains the strongest of these:
- 2.4979888214

but does not retain the other two inside its predicted Top-42.

Therefore normalized regret remains exactly equal to B2 even though discrete
Top-1% overlap falls from 3 to 1.

## Preregistered decision checks

PASS:
- full-test Spearman beats B2
- clean-region Spearman beats B2
- normalized regret not worse than B2
- full-test Spearman beats clean pair expansion
- positive Spearman at mutation orders 3, 4 and 5
- residual-training OOF Spearman positive
- deterministic rank replay

FAIL:
- Top-1% hits >= B2

Observed:
- 7 / 8 checks pass

Strict preregistered decision:

    V9_CONTEXTUAL_ESM_RESIDUAL_REJECT

The gate is not changed post hoc.

## Architectural conclusion

The contextual representation is materially better than context-free pair
composition and is the first V9 development mechanism to exceed B2 in global
rank correlation and in the directly affected clean-supported region.

However, because the preregistered elite-overlap gate failed, this exact
architecture is not retained as the final V9 architecture.

Do not add scalar calibration, mutation-count scaling, clipping, or a post-hoc
Top-K objective to rescue this run.

The next genuine architectural source, per protocol, is structure-aware
relational information while preserving the successful contextual sequence
signal as a diagnostic reference.

IRED remains development-only. Phase 3 remains closed.
