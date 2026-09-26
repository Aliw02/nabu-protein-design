# B2-Residual Context Development Protocol

Status: POST-PHASE-2 DEVELOPMENT
Dataset: FLIP2 TrpB two-to-many (already revealed development evidence)
Phase 3: CLOSED
Scientific claim: false

## Motivation

The frozen TrpB blind result showed that the previous clean-singleton additive
baseline was badly mis-scaled on higher-order variants. The contextual model
therefore spent much of its capacity cancelling that baseline rather than
adding independent information beyond B2.

This experiment changes exactly one conceptual rule:

    OLD:
      clean-additive baseline + Context(clean-additive residual)

    NEW:
      B2 baseline + Context(B2 residual)

No structure channel, no Pareto consensus, no mutation-count weighting, no
hand-coded residue classes, no clipping, no test-time tuning.

## Training target

Fit the frozen NABU V8.3 model on the fit/train partition only.

Use its cross-fitted B2 prediction already produced by the model:

    b2_oof[i]

Define:

    residual_target[i] = y[i] - b2_oof[i]

The frozen ESM-2 contextual feature is unchanged:
- facebook/esm2_t6_8M_UR50D
- revision c731040fcd8d73dceaa04b0a8e6329b345b0f5df
- 640 dimensions:
  - global mean final-hidden delta vs WT (320)
  - mean mutation-site final-hidden delta vs WT (320)

The readout is unchanged:
- StandardScaler
- RidgeCV
- alphas: 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1000, 10000
- deterministic 5-fold CV
- no validation/test targets used for fit

A deterministic shuffled-residual control uses the same features and readout
with residual targets permuted using seed 161.

## Higher-order scoring

For every test sequence:

    B2_MAIN_ONLY = frozen B2 score

    B2_PLUS_CONTEXT_RESIDUAL =
        B2_MAIN_ONLY + ContextResidual(sequence)

    B2_PLUS_SHUFFLED_CONTEXT =
        B2_MAIN_ONLY + ShuffledContextResidual(sequence)

No clean-singleton coverage gate is used.

## Development acceptance gate

The architecture is eligible to freeze for a new untouched benchmark only if
ALL are true on the already-revealed TrpB development test:

1. full-test Spearman(B2 + Context) > Spearman(B2)
2. Top-1% hits(B2 + Context) >= Top-1% hits(B2)
3. normalized Top-1% regret(B2 + Context) <= B2
4. Spearman is positive for mutation order 3
5. Spearman is positive for mutation order 4
6. Spearman(B2 + Context) > Spearman(B2 + shuffled Context)
7. deterministic replay/hash contracts pass

If this gate fails, do not consume a new untouched benchmark.

If it passes, freeze this exact architecture before accessing any new
untouched benchmark target labels.
