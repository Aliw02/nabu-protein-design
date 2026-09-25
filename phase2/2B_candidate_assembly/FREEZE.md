# Phase 2B Freeze

## Status

**PASS WITH LIMITATION**

Candidate Assembly V1 is operational, reference-aware, deterministic and leak-resistant.

## Engineering result

RhlA 5% development POC:

- reference sequence inferred consistently from mutation identities: `SFQRAQLSQHA`;
- validated mutation vocabulary: 30 substitutions;
- 64 triple-mutant proposals;
- 64 / 64 proposals scoreable;
- 0 measured-overlap proposals;
- 0 duplicate proposals;
- hidden/unmeasured truth perturbation leaves proposals unchanged;
- proposal list frozen and hashed before virtual truth lookup.

Of the 64 assembled proposals:

- 26 existed in the historical RhlA table;
- 38 were novel-to-table combinations with unknown virtual truth.

This proves assembly is not simply re-ranking the supplied candidate table.

## Initial quality limitation

At 5% measured evidence, the historical matched subset was weak:

- best known fitness = -0.228583;
- Top-1% hits = 0.

Therefore engineering validity alone is not sufficient to enable assembly early in a campaign.

## Evidence-maturity ablation

The assembler itself was held fixed.

Only measured evidence changed:

| Evidence | Pair entries | Triple frontier scoreability | Matched best | Top-1% hits | Mean percentile |
| --- | ---: | ---: | ---: | ---: | ---: |
| 5% | 169 | 2299 / 2757 | -0.228583 | 0 | 0.746747 |
| 10% | 250 | 2605 / 2690 | -0.228583 | 0 | 0.667226 |
| 20% | 333 | 2537 / 2537 | 0.308292 | 2 | 0.792204 |

The improvement is not monotonic at low evidence.

At the 20% RhlA state, assembly quality improved materially.

## Important architecture finding

The improvement at 20% occurred while:

- router mode remained `B3_PROTECTED_NO_HIGHER_ORDER`;
- triplet memory entries = 0;
- quartet memory entries = 0.

Therefore the observed assembly improvement came from a denser pairwise evidence base, not from higher-order memory activation.

## 2C decision

2C may begin, but assembly must not be enabled blindly from the first round.

Carry forward an explicit **maturity-gated assembly** principle.

Do not freeze 20% as a universal biological threshold.

The 20% RhlA result is development evidence showing that assembly usefulness depends on evidence maturity.

2C should test a generic evidence-state gate rather than a dataset-name or fixed-percentage rule.
