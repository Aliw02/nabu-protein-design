# NABU V9 Contextual Interaction — Frozen ESM-2 Residual Model

Status: POST-PHASE-2 DEVELOPMENT
Phase 3: CLOSED
Scientific claim: false

## Purpose

Replace the rejected context-free higher-order pair sum with one global
context-dependent interaction correction.

The previous clean-pair diagnostics established:

- clean pair residual is learnable on unseen pair edges;
- affine calibration repairs pair-residual magnitude;
- summing calibrated pair effects in higher-order variants is not sufficient;
- the summed pair signal still correlates with the true missing higher-order
  residual, so interaction information exists but composition is context dependent.

This experiment changes the interaction representation, not Phase 2.

## Frozen external encoder

Use the public pretrained ESM-2 checkpoint:

    facebook/esm2_t6_8M_UR50D

Pinned revision:

    c731040fcd8d73dceaa04b0a8e6329b345b0f5df

The encoder is frozen. No IRED target is used to update ESM-2 weights.

## Training evidence

Use only the IRED official fit/train partition:

    set == train AND validation != True

No validation or test target is used for fitting the contextual model.

Construct clean additive evidence:

- WT target;
- 1,037 measured singleton effects;
- double-mutant rows only when both singleton counterparts exist.

For WT and measured singles, contextual interaction residual is defined as 0.

For clean doubles:

    residual =
        y(double)
        - y(single_i)
        - y(single_j)
        + y_WT

Thus the contextual residual training set contains:

- WT: residual 0
- all measured singles: residual 0
- all clean doubles: measured clean pair residual

The 949 doubles lacking one or both singleton counterparts are excluded from
contextual residual fitting and reported.

## Context representation

For every sequence, obtain the frozen final-layer ESM-2 hidden state.

Let R be the WT/reference sequence embedding and S the candidate embedding.

Two 320-dimensional context channels are used:

1. GLOBAL_DELTA
   Mean residue embedding of S minus mean residue embedding of R.

2. MUTATION_SITE_DELTA
   Mean over mutated positions of:
       hidden_S[position] - hidden_R[position]

   For WT, this channel is the all-zero vector.

Final feature vector:

    concat(GLOBAL_DELTA, MUTATION_SITE_DELTA)

Dimension: 640.

This representation depends on the complete candidate sequence. Therefore the
same mutation can receive a different contextual representation when other
mutations are present.

No amino-acid class lists, motif rules, hand-authored biochemical properties,
or target-specific positional rules are used.

## Residual readout

Fit one linear Ridge readout from the 640-dimensional frozen contextual feature
to the interaction residual.

Preprocessing:
- StandardScaler fitted on the fit/train residual-training rows only.

Regularization is selected internally from the fixed generic grid:

    1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1000, 10000

using deterministic 5-fold cross-validation within the residual-training set.

No IRED validation/test labels participate in alpha selection.

## Higher-order scoring

For a test candidate where every mutation has a measured singleton:

    score =
        y_WT
        + sum(clean_main(m))
        + contextual_residual(sequence)

If any mutation lacks a measured singleton, return the unchanged B2 score for
that entire candidate and mark it as B2 fallback.

There is:
- no pair summation;
- no mutation-count scaling;
- no clipping;
- no post-hoc calibration on validation/test;
- no threshold selection.

## Comparison arms

1. PHASE2_V8_3_ABSTENTION
2. B2_MAIN_ONLY
3. V9_1_DEGREE_NORMALIZED_TRANSFER
4. V9_CLEAN_CALIBRATED_PAIR_EXPANSION
5. V9_CONTEXTUAL_ESM_RESIDUAL

IRED test is development-only because its targets were already revealed.

## Primary metrics

- full-test Spearman;
- NDCG;
- Top-1% hits / recall / enrichment;
- normalized Top-1% regret;
- complete-clean-main-region metrics;
- metrics by mutation count;
- selected Ridge alpha;
- residual-training OOF Spearman / RMSE;
- feature/prediction hashes;
- deterministic rank replay.

## Decisive rule

Retain V9_CONTEXTUAL_ESM_RESIDUAL as the V9 development architecture only if ALL:

1. full-test Spearman > B2_MAIN_ONLY;
2. complete-clean-main-region Spearman > B2 on the same rows;
3. Top-1% hits >= B2_MAIN_ONLY;
4. normalized Top-1% regret <= B2_MAIN_ONLY;
5. full-test Spearman > V9_CLEAN_CALIBRATED_PAIR_EXPANSION;
6. Spearman > 0 separately for mutation orders 3, 4, and 5;
7. residual-training OOF Spearman > 0;
8. deterministic replay preserves identical rank ordering.

If this fails, do not create scalar/normalization patch versions. The next
architecture must add genuinely richer relation information, with structure
being the next candidate source.

## Boundary

This experiment is development evidence only.
A retained model still requires a new untouched external benchmark before any
Phase-2-style scientific claim can be revisited.
Phase 3 remains closed.
