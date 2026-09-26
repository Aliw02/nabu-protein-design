# NABU V9 Simple-Architecture Ablation

Status: POST-PHASE-2 DEVELOPMENT
Phase 3: CLOSED
Scientific claim: false

## Question

Can the remaining elite-ranking failure be solved by either:

A. preserving two complementary rankings instead of collapsing them to one fitness scalar; or
B. adding one simple contact-aware relational channel to the successful frozen-ESM contextual residual model?

No scalar rescue, mutation-count scaling, clipping, target-specific rules, or post-hoc Top-K tuning is allowed.

## Shared data boundary

Use the pinned IRED source and official fit/train partition exactly as in prior V9 development.

- fit/train targets may be used for fitting;
- validation targets are not used for fitting or model selection;
- test targets are development evaluation only because Phase 2D already revealed them.

Frozen V8.3 remains unchanged.
Phase 3 remains closed.

## Arm A — B2 + Context Consensus

Inputs:
- B2_MAIN_ONLY score
- previously defined V9_CONTEXTUAL_ESM_RESIDUAL score

No training is added.

For all test candidates:

1. Convert both scores independently to descending fractional ranks.
2. Perform deterministic two-objective non-dominated sorting:
   a candidate dominates another only when it is at least as good in both ranks
   and strictly better in at least one.
3. Assign Pareto layers from best to worst.
4. Within each layer, tie-break by the arithmetic mean of the two fractional ranks.
5. Convert the final deterministic ordering to a monotone score for metric evaluation.

There are no weights and no fitted parameters.

## Arm B — Contact-aware contextual residual

Use the same frozen encoder and revision as the contextual V9 experiment:

    facebook/esm2_t6_8M_UR50D
    revision c731040fcd8d73dceaa04b0a8e6329b345b0f5df

The encoder remains frozen.

Retain the existing 640-dimensional sequence-context feature:

- global mean final-hidden delta vs WT: 320 dims
- mutation-site mean final-hidden delta vs WT: 320 dims

Add one 320-dimensional contact-neighborhood channel.

For candidate residue hidden deltas D[q] and the frozen ESM predicted contact
probability matrix C[p,q]:

For every mutated position p:

    neighborhood(p) =
        sum_q C[p,q] * D[q] / max(sum_q C[p,q], epsilon)

Then:

    CONTACT_NEIGHBOR_DELTA =
        mean_p neighborhood(p)

For WT, this channel is the all-zero vector.

Final feature dimension:

    320 + 320 + 320 = 960

This uses no residue-class lists, motif rules, manually chosen positions,
distance cutoffs, or target-specific biochemical features.

Fit exactly the same RidgeCV residual readout and alpha grid used by the
sequence-context model:

    1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1000, 10000

with deterministic 5-fold fit/train-only cross-validation.

Higher-order score on complete-clean-main rows:

    clean_additive_baseline + predicted_contact_context_residual

Rows missing any measured singleton use the unchanged B2 score for the entire candidate.

## Comparison

Evaluate:

1. B2_MAIN_ONLY
2. V9_CONTEXTUAL_ESM_RESIDUAL
3. V9_CONSENSUS_B2_CONTEXT
4. V9_CONTACT_AWARE_CONTEXT

Report:
- full-test Spearman
- NDCG
- Top-1% hits / recall / enrichment
- normalized Top-1% regret
- complete-clean-main-region metrics
- mutation-order metrics
- residual-training OOF metrics for the contact-aware arm
- deterministic hashes

## Arm acceptance

An arm is individually acceptable only if ALL are true:

1. full-test Spearman >= B2_MAIN_ONLY;
2. complete-clean-main Spearman >= B2 on the same rows;
3. Top-1% hits >= B2_MAIN_ONLY;
4. normalized Top-1% regret <= B2_MAIN_ONLY;
5. positive Spearman at mutation orders 3, 4, and 5;
6. deterministic replay is identical.

For Contact-aware only:
7. residual-training OOF Spearman > 0.

## Predeclared ablation winner

- If exactly one arm passes its acceptance gate, that arm wins.
- If both pass, choose the arm with higher full-test Spearman.
- If full-test Spearman ties to 1e-12, choose higher complete-clean-main Spearman.
- If still tied, choose higher Top-1% hits.
- If neither passes, the ablation result is NO_SIMPLE_SOLUTION_PASS.

The gate and winner rule are not changed after observing test metrics.
