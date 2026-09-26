# NABU V9 TrpB Two-to-Many Blind External Validation

Status: PREREGISTERED / BLIND EXTERNAL VALIDATION
Phase 3: CLOSED
Scientific claim before reveal: false

## Frozen architecture

The architecture under test is the IRED-development winner:

    V9_CONSENSUS_B2_CONTEXT

Frozen components:

1. V8.3 B2 main-only evidence score.
2. Clean additive baseline from WT + measured singleton evidence.
3. Frozen ESM-2 contextual residual:
   - facebook/esm2_t6_8M_UR50D
   - revision c731040fcd8d73dceaa04b0a8e6329b345b0f5df
   - 640-dimensional context feature:
       a. global mean final-hidden delta vs WT (320)
       b. mean mutation-site final-hidden delta vs WT (320)
   - RidgeCV alpha grid:
       1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1000, 10000
   - deterministic 5-fold train-only model selection.
4. Parameter-free two-objective Pareto consensus over:
       B2 ranking
       Context ranking
   with arithmetic mean fractional rank only as a within-Pareto-layer tie-break.

No TrpB test target can modify the architecture, feature definition, model
selection, fallback rule, consensus rule, thresholds, or pass gate.

Canonical IRED development evidence:
- V9 simple ablation canonical run: 36238279578
- architecture commit used by that run: eeabc2ef4811abf64776fff06eda84925aafb3ba
- recorded result summary commits:
  - f6d6ad1375ac299fc48b99364cf8677354d0b5c9
  - f68ae43a5d55a0edef88950298c875d7debecc07

## External benchmark

FLIP2 TrpB two-to-many, versioned Zenodo record:

    https://zenodo.org/records/18433203

File:

    trpb/two_to_many.csv.gz

Published file MD5:

    a611408d2db06907a023ddc8a1d94c12

Expected official split counts:
- fit/train: 8,633
- validation: 2,158
- test: 217,507

Published split semantics:
- train: WT, single mutants, and double mutants;
- test: higher-order variants.

The workflow computes and freezes the SHA256 of the exact downloaded source.

## Blind isolation contract

The workflow uses separate artifacts.

### Prediction-visible artifact

Contains:
- fit/train rows WITH targets;
- test rows WITHOUT target;
- source/split manifest.

It does NOT contain test targets.

### Reveal-only artifact

Contains:
- test row_id;
- test target.

Prediction and training jobs never download the reveal artifact.

The sequence is:

    source download
      -> split and mask
      -> train frozen model on fit/train only
      -> shard test prediction
      -> aggregate predictions
      -> freeze prediction file + hashes
      -> ONLY THEN download reveal artifact
      -> evaluation

Validation targets are not used for fitting or model selection.

## Deterministic duplicate-evidence adapter

TrpB combines multiple sub-landscapes and may contain repeated evidence for the
same mutation identity. Before seeing any test target, the following generic
adapter is frozen:

- V8.3 B2 receives all fit/train rows unchanged.
- WT target for clean-additive context is the arithmetic mean over fit/train WT rows.
- Singleton target for each mutation identity is the arithmetic mean over all
  fit/train singleton rows carrying that identity.
- Each eligible double-mutant training row retains its own measured target and
  uses the frozen WT/singleton means to define its clean interaction residual.
- No test or validation target participates in this adapter.

This adapter has no TrpB-specific amino-acid, position, motif, or fitness rule.

## Contextual residual training target

For WT and singleton fit/train rows:

    residual = 0

For every double-mutant fit/train row whose two singleton identities are both
available:

    residual =
        y(double)
        - mean_y(single_i)
        - mean_y(single_j)
        + mean_y(WT)

Double rows lacking one or both singleton identities are excluded from the
contextual residual fit and counted.

## Test scoring

For a test variant:

B2:
- unchanged V8.3 B2 semantics.

Context:
- if all mutation identities have clean singleton evidence:

      clean_additive =
          mean_y(WT) + sum(clean_main(m))

      context_score =
          clean_additive + frozen_contextual_residual(sequence)

- otherwise:

      context_score = B2 score

Consensus:
- rank all 217,507 B2 scores;
- rank all 217,507 Context scores;
- compute exact two-objective Pareto layers;
- within a layer, tie-break by mean fractional rank;
- candidate row_id is the final deterministic tie-break.

The large-N implementation is an O(N log N) exact implementation of the same
two-objective Pareto semantics used on IRED. Unit tests require exact equality
with the original O(N^2) implementation on controlled data.

## Frozen predictions before reveal

Before reveal the workflow writes:

- TRPB_FROZEN_PREDICTIONS.csv.gz
- TRPB_FREEZE_MANIFEST.json

The manifest contains:
- source SHA256;
- row count;
- semantic prediction SHA256;
- B2 rank SHA256;
- Context rank SHA256;
- Consensus rank SHA256.

The reveal/evaluation job runs only after this artifact exists.

## Evaluation

Report for:
- B2_MAIN_ONLY
- V9_CONTEXTUAL_ESM_RESIDUAL
- V9_CONSENSUS_B2_CONTEXT

Metrics:
- Spearman
- NDCG
- Top-1% hits / recall / enrichment
- normalized Top-1% regret
- best true target in predicted Top-1%
- complete-clean-main region metrics
- mutation-order metrics
- coverage
- frozen hashes

## Published FLIP2 reference points

For TrpB two-to-many, published FLIP2 reference values include:

- Ridge one-hot:
    Spearman 0.437
    NDCG 0.994

- CARP-640M supervised:
    Spearman 0.509
    NDCG 0.996

These values are comparison references only and are not used for fitting.

## Preregistered decisions

### Internal external-transfer gate

PASS only if all are true:

1. Consensus full-test Spearman > B2 full-test Spearman.
2. Consensus complete-clean-main Spearman > B2 on the same rows.
3. Consensus Top-1% hits >= B2 Top-1% hits.
4. Consensus normalized Top-1% regret <= B2.
5. Consensus Spearman is positive for mutation order 3.
6. Consensus Spearman is positive for mutation order 4.
7. Frozen prediction hashes are present and deterministic.
8. Official source/split contracts pass.
9. Validation targets are unused for fitting.

### Published-benchmark gate

PASS only if:
- the internal external-transfer gate passes; and
- Consensus full-test Spearman >= 0.437 (published Ridge one-hot reference).

### Strong published benchmark

Reported, not required for the primary gate:

    Consensus Spearman >= 0.509

## Interpretation

If both the internal external-transfer gate and the published-benchmark gate
pass, this is the first new untouched external evidence that the V9 consensus
architecture repairs the higher-order generalization failure exposed in Phase 2D.

It still does not retroactively convert the old Phase-2D IRED result into a pass.

If the internal gate fails, do not tune on TrpB test targets.
If the internal gate passes but the published-benchmark gate fails, report
transfer success but insufficient benchmark competitiveness.

Phase 3 remains closed until the external result is reviewed and explicitly frozen.
