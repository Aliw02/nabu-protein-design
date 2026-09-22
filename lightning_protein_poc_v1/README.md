# Lightning Protein/DMS POC V1

Input: 200 KCNH2 single-mutant DMS rows across positions 545-555.

Blind protocol:
- Deterministic 5-fold holdout within each position (seed 161).
- Hidden labels are excluded from prediction.
- Each fold's predictions are serialized and SHA-256 frozen before truth reveal.

Arms:
- B0_POSITION_MEAN: visible mean for the position.
- B1_ALT_GLOBAL: visible global mean for the target amino acid.
- CORE_ADDITIVE: position mean + target-amino-acid residual.
- LIGHTNING_DREAM: position mean + relative-resonance-weighted residual from the same target amino acid at other visible positions.

Lightning-inspired concepts:
- Typed contract: CORE / IDENTITY_ANCHOR / QUALIFIER / QUANTITY / CONTEXT.
- Exact identity only; no fuzzy aliasing.
- Relative resonance learned from visible position response profiles.
- Dream candidates are hidden/unmeasured variants; they are retained as WARM or DORMANT.
- Freeze-before-reveal via SHA-256.
- No amino-acid property table or manual protein-domain vocabulary.

This is an engineering POC, not evidence of biological validity or laboratory performance.
