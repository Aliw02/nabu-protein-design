# NABU B2-Residual Context — NucB Two-to-Many Untouched Blind Protocol

Status: PREREGISTERED / NOT YET RUN
Dataset: FLIP2 NucB two-to-many
Phase 3: CLOSED
Scientific claim before reveal: false

## Untouched status

Repository search before protocol creation found no prior NucB experiment,
adapter, prediction, result, or target use in this project.

The benchmark file is pinned to the official FLIP2 Zenodo record:

    https://zenodo.org/records/18433203

File:

    nucb/two_to_many.csv.gz

Published MD5:

    cdedaadd57f2b148950b357e9f24c238

NucB measures nuclease activity at pH 7 and FLIP2 bins activity into four
ordered activity levels. The official two-to-many split trains on variants
with 0, 1, or 2 mutations and tests on higher-order variants.

## Architecture freeze

This blind run is permitted only if the preregistered TrpB development gate
for B2_PLUS_CONTEXT_RESIDUAL passes.

If permitted, the architecture is frozen exactly as:

    B2_MAIN_ONLY
    +
    frozen ESM-2 contextual prediction of B2 cross-fit residual

Training target:

    residual_target[i] = y[i] - b2_oof[i]

Context encoder:
- facebook/esm2_t6_8M_UR50D
- revision c731040fcd8d73dceaa04b0a8e6329b345b0f5df
- frozen
- 640-dimensional feature:
  - global mean final-hidden delta vs WT (320)
  - mean mutation-site final-hidden delta vs WT (320)

Readout:
- StandardScaler
- RidgeCV
- alpha grid:
  1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1000, 10000
- deterministic 5-fold CV
- no validation/test target used for fit

Higher-order score:

    score = B2(sequence) + ContextResidual(sequence)

No:
- Pareto consensus
- clean-singleton additive baseline
- mutation-count scaling
- target clipping
- residue wordlists
- structural thresholds
- NucB-specific tuning
- validation-target tuning
- post-reveal retuning

A deterministic shuffled-residual control is retained as a negative control,
but it cannot influence model selection.

## Blind isolation

The prepare job may access the source file only to:
1. verify the exact MD5;
2. assign fit/validation/test rows from the official split columns;
3. remove test targets from the prediction-visible artifact;
4. seal test row_id + target into a reveal-only artifact;
5. record split and mutation-order metadata.

Prediction/training jobs cannot download the reveal artifact.

Sequence/reference scope checks occur before training:
- fixed sequence length is required;
- fit/train mutation orders must be <=2;
- test mutation orders must be >=3.

If these contracts fail, the result is ADAPTER_SCOPE_FAIL, not a model failure.

Predictions are frozen and hashed before the reveal job can download targets.

## Metrics

Because NucB targets are four ordered activity bins, the primary ranking
metrics are:
- Spearman
- NDCG

Also report:
- mutation-order Spearman/NDCG
- Top-1% metrics as descriptive only (ties make them less suitable as a gate)
- B2 support coverage
- prediction hashes
- deterministic replay contracts

Arms:
1. B2_MAIN_ONLY
2. B2_PLUS_CONTEXT_RESIDUAL
3. B2_PLUS_SHUFFLED_CONTEXT

## Published references

FLIP2 NucB two-to-many published references:
- Ridge one-hot: Spearman 0.606, NDCG 0.943
- Ridge one-hot + likelihoods: Spearman 0.623, NDCG 0.946
- CARP-640M supervised: Spearman 0.717, NDCG 0.970
- ESMC-300M supervised: Spearman 0.723, NDCG 0.966

The primary published benchmark threshold is the plain Ridge one-hot
Spearman 0.606. The 0.723 ESMC result is reported as a strong reference,
not a required primary gate.

## Preregistered blind PASS gate

PASS requires ALL:

1. Spearman(B2 + Context) > Spearman(B2)
2. NDCG(B2 + Context) >= NDCG(B2)
3. Spearman(B2 + Context) > Spearman(B2 + shuffled Context)
4. Spearman(B2 + Context) >= 0.606
5. full prediction freeze completed before reveal
6. source MD5 and split semantics pass
7. validation targets unused for fit
8. all predictions finite and deterministic
9. Spearman is positive for every test mutation order with at least 100 rows

Strong-reference result is separately reported if:

    Spearman(B2 + Context) >= 0.723

If the primary gate passes, this is new untouched external evidence that the
post-Phase-2 B2-residual Context architecture transfers beyond IRED and TrpB.

No NucB target may be used for retuning after reveal.
