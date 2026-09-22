# Lightning GB1 Blind Combinatorial Design V3

## Result

Primary candidate lists were frozen before hidden candidate fitness was revealed.

| Arm | Top-1 | Mean Top-5 | Mean Top-10 | Median Top-10 | > WT (1.0) | Best recovered |
|---|---:|---:|---:|---:|---:|---:|
| Lightning Resonance | 3.4540 | 3.7090 | 3.5541 | 3.4267 | 90% | 6.3435 |
| Ridge main+pair | 4.4532 | 3.8931 | 3.6407 | 3.2491 | 100% | 6.3435 |
| Additive Memory | 2.5417 | 2.7397 | 2.1676 | 1.8528 | 80% | 4.4532 |
| Random hash | 0.0090 | 0.0036 | 0.0370 | 0.0034 | 0% | 0.3367 |
| Shuffled-label Lightning | 0.0049 | 0.0060 | 0.0103 | 0.0058 | 0% | 0.0450 |

Lightning mean Top-10 is **1.64x additive**, and the deterministic shuffled-label control falls by **99.7%** relative to true-label Lightning. Lightning recovered `WWFG` with measured fitness **6.3435** at rank 5.

The strongest learned baseline, Ridge main+pair, remains slightly higher on mean Top-10: Lightning is **97.6%** of Ridge. Therefore this is strong blind compositional-design evidence, but not a claim that Lightning dominates a strong interaction regression baseline.
