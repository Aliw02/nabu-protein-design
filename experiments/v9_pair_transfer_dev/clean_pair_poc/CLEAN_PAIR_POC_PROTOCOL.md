# Clean Main → Clean Pair Residual → Unseen-Pair Generalization POC

Status: POST-PHASE-2 DEVELOPMENT DIAGNOSTIC
Phase 3: CLOSED
Purpose: determine whether the V9 failure comes from contaminated main effects or from insufficient pair representation.

## Question

Can pair epistasis be predicted for an unseen double-mutant pair after defining
main effects only from measured single mutants relative to WT?

This POC does not modify V8.3, V9.0, V9.1, Phase 2, or any frozen scientific claim.

## Dataset role

FLIP2 IRED two-to-many is now development-only because its targets were already
revealed during Phase 2D.

Only the official fit/train partition is used for model fitting and cross-validation:

    set == train AND validation != True

The official validation and test partitions are not used by this POC.

Expected fit composition:
- 1 WT
- 1,037 single mutants
- 2,708 double mutants

## Clean main effects

Let y_WT be the measured WT target.

For a mutation m with a measured single-mutant row:

    main(m) = y(single_m) - y_WT

No double-mutant or higher-order row contributes to main(m).

## Clean pair residual

For a double mutant (i,j) where both single mutations are measured:

    epsilon(i,j) = y(i,j) - [y_WT + main(i) + main(j)]

Equivalently:

    epsilon(i,j) = y(i,j) - y(single_i) - y(single_j) + y_WT

Double rows missing either singleton counterpart are structurally excluded from
the CLEAN_PAIR subset and reported separately. This exclusion is identity-only,
not target-aware.

## Deterministic unseen-pair cross-validation

Each clean pair edge is assigned to one of five folds by:

    fnv1a32(candidate_id + "|NABU_CLEAN_PAIR_CV|161") % 5

For each fold:
- train on the other four folds;
- predict clean pair residuals for held-out pair edges;
- no held-out residual enters fitting.

Held-out edges are categorized as:
1. NODE_COVERED: both mutation identities occur in at least one training edge;
2. ONE_NODE_COLD;
3. TWO_NODE_COLD.

Primary generalization analysis uses NODE_COVERED unseen edges.
Cold-node cases are reported separately.

## Models / controls

A. ZERO_INTERACTION
    epsilon_hat = 0

B. TRAIN_MEAN
    epsilon_hat = mean(train epsilon)

C. NODE_ADDITIVE_RELATION
    epsilon_hat(i,j) = bias + u_i + u_j

The node-additive relation is fitted by deterministic minimum-norm least squares
on training clean-pair residuals only. No regularization strength, rank, threshold,
or post-hoc hyperparameter is tuned.

D. PERMUTED_NODE_ADDITIVE_CONTROL
    Same model and graph, but training residuals are deterministically permuted
    with seed 161 inside each fold.

## Metrics

For NODE_COVERED unseen edges:
- Spearman correlation between epsilon_hat and true epsilon;
- Pearson correlation;
- RMSE;
- MAE;
- sign accuracy;
- improvement in RMSE over ZERO_INTERACTION.

Also report:
- clean-pair coverage;
- node coverage by fold;
- train design rank;
- residual distribution;
- fold-level metrics;
- pooled out-of-fold metrics;
- deterministic prediction hashes.

## Decision

The current identity-only pair-transfer hypothesis is considered supported only if
NODE_ADDITIVE_RELATION, on NODE_COVERED unseen edges:

1. has positive Spearman;
2. beats ZERO_INTERACTION on RMSE;
3. beats TRAIN_MEAN on RMSE;
4. beats the permuted control on Spearman;
5. sign accuracy > 0.50;
6. the direction is consistent in at least 4 of 5 folds for Spearman improvement
   over the permuted control.

If this fails, do not continue scalar/normalization tuning of V9 pair transfer.
The next architecture must use a richer relational representation such as
sequence-context or structure-aware pair features.

## Boundary

A positive result is development evidence only.
A negative result is a diagnostic limitation, not a new Phase-2 result.
