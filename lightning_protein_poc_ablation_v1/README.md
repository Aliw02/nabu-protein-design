# Lightning Protein/DMS POC — Ablation V1

Dataset: 200 KCNH2 single-mutant DMS variants, positions 545–555.

## Protocol

- Reuses the exact frozen 5-fold mutant assignment from POC V1.
- Hidden predictor input contains descriptors only; hidden DMS truth is absent.
- Each fold prediction payload is SHA-256 frozen before truth merge.
- No amino-acid property table, protein-domain vocabulary, embeddings, or hidden-label tuning.
- Relative resonance is learned only from visible response profiles.

## Arms

- **B0_GLOBAL_MEAN** — visible global mean.
- **B1_A4_CORE_POSITION** — typed CORE/context represented by the visible mean at the exact position.
- **B2_A4_PLUS_MEMORY** — B1 plus position-centered ALT residual memory learned across visible positions.
- **B3_PLUS_RELATIVE_RESONANCE** — B2 upgraded by mutual row-relative resonance transfer when visible same-ALT support exists; B2 fallback otherwise.
- **B4_DREAM_WARM_ONLY** — same score as B3 but only WARM candidates are retained.
- **B5_FULL_DREAM_DORMANT** — WARM candidates plus unsupported candidates retained DORMANT with the B2 estimate.

## Critical result

B2 Memory Spearman: **0.463645**  
B5 Full Spearman: **0.525015**  
Delta: **+0.061370**

B2 ROC-AUC: **0.736000**  
B5 ROC-AUC: **0.754900**  
Delta: **+0.018900**

B2 mean Top-3 hidden score: **59.3647**  
B5 mean Top-3 hidden score: **59.6448**

Per-position Spearman sign-flip p-value (B5 vs B2): **0.176758**.  
Cluster bootstrap 95% CI for Spearman delta: **[-0.038928, 0.177809]**.

Dream states: **197 WARM / 3 DORMANT**.

## Interpretation boundary

This is an engineering ablation on one small DMS slice. It tests whether the Lightning-inspired representation/memory/resonance mechanics add blind predictive ranking signal. It does **not** establish biological validity, de-novo protein design performance, or laboratory efficacy.
