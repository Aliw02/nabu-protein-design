# Lightning Protein Design POC — Cross-Assay Replication V2

## Decision

**Engineering transfer POC: PASS. Protein-design capability: not yet proven.**

Frozen Lightning-style Memory + Relative Resonance improved pooled out-of-fold ranking over the Memory-only arm in **5/5 assays**. Deterministic visible-label shuffle controls collapsed to approximately zero in **5/5 assays**, and whole-position holdout memory transfer remained positive in **5/5 assays**. However, top-1 candidate selection is not uniformly improved, so this is evidence of transferable ranking/memory structure, not yet evidence that the architecture can reliably design the best unseen protein variant.

## Cross-assay results

| Dataset | Rows | B2 Spearman | B5 Spearman | Δ B5-B2 | B5 group ρ | ALT group ρ | B5 Top-1 pct | ALT Top-1 pct | Shuffle B5 ρ | Position-holdout Memory ρ |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| KCNH2_subset | 200 | 0.464 | 0.525 | +0.061 | 0.322 | 0.338 | 0.726 | 0.750 | -0.008 | 0.276 |
| GFP_AEQVI_Sarkisyan_2016 | 1084 | 0.481 | 0.505 | +0.025 | 0.185 | 0.138 | 0.792 | 0.748 | -0.045 | 0.131 |
| SUMO1_HUMAN_Weile_2017 | 1700 | 0.585 | 0.592 | +0.007 | 0.128 | 0.113 | 0.695 | 0.704 | -0.001 | 0.084 |
| UBC9_HUMAN_Weile_2017 | 2563 | 0.672 | 0.690 | +0.018 | 0.191 | 0.162 | 0.714 | 0.707 | -0.007 | 0.149 |
| PABP_YEAST_Melamed_2013 | 1187 | 0.748 | 0.774 | +0.025 | 0.376 | 0.291 | 0.812 | 0.770 | 0.017 | 0.204 |

## Aggregate

- Total evaluated variants: **6,734** across **5 assays**.
- B5 pooled Spearman > B2: **5/5 assays**.
- Mean pooled Spearman gain from Relative Resonance/Dream over Memory: **+0.0271**.
- Shuffle controls |ρ| < 0.05: **5/5 assays**.
- Whole-position holdout Memory transfer positive: **5/5 assays**.
- B5 mean fold-position group Spearman > ALT-global baseline: **4/5 assays**.
- B5 Top-1 percentile > ALT-global baseline: **3/5 assays**.

## Paired design-selection diagnostics

| Dataset | B5-B2 group-ρ sign p | B5-B2 Top-1 truth Δ | Top-1 sign p |
|---|---:|---:|---:|
| KCNH2_subset | 1 | -1.297 | 0.4807 |
| GFP_AEQVI_Sarkisyan_2016 | 1 | -0.002202 | 0.06402 |
| SUMO1_HUMAN_Weile_2017 | 0.567 | +0.009585 | 0.7338 |
| UBC9_HUMAN_Weile_2017 | 0.04241 | +0.005817 | 0.514 |
| PABP_YEAST_Melamed_2013 | 0.007178 | +0.02939 | 7.011e-05 |

Strongest fresh result is PABP: B5 improved pooled Spearman from 0.748 to 0.774; fold-position group Spearman improved from the B2 arm, with sign-test p≈0.0072, and Top-1 selected hidden truth improved with p≈7.0e-5. UBC9 also shows a group-Spearman advantage (p≈0.042). GFP and SUMO1 show smaller/non-significant B5-vs-B2 group gains, and KCNH2 does not show a Top-1 gain.

## What the controls say

The visible-label shuffle is the most important sanity check here: after breaking the mapping between visible mutations and visible measurements, B5 pooled Spearman is approximately zero on every assay. That argues against the observed signal being caused by the mutation IDs, fold assignment, or static sequence descriptors alone.

The whole-position holdout test removes every measurement from a target position. The Memory arm still has positive Spearman on all five assays, indicating that the ALT residual memory transfers some information across positions. Relative Resonance cannot fully operate for a completely unseen position because that position has no visible response profile; this stress test therefore isolates the cross-position Memory component.

## Boundary

This does **not** establish de novo protein design. These experiments rank hidden measured single substitutions. The next scientific step is a new, prospectively frozen candidate-generation test where Lightning proposes mutations not used in the benchmark ranking experiment and an independent evaluator scores them. The current revealed assays must now be treated as development data, not reused as a fresh confirmation set.