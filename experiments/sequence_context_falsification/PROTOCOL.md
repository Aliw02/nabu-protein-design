# NABU Sequence-Context Falsification 01

Date: 2026-10-01
Branch: exp/nabu-sequence-context-falsification
Status: preregistered discriminative test
NucB: FORBIDDEN
Phase 3: CLOSED
Governance: SCIENTIFIC_CONSTRAINTS.md on main

## Scientific question

Does a frozen pretrained protein language model provide candidate-conditioned sequence-context information about higher-order fitness residuals that is not already explained by low-order O1/O2/interaction-shape statistics or intrinsic single-mutation plausibility?

This is a falsification experiment, NOT an architecture-development experiment.

## Competing hypotheses

### H_SEQ_CONTEXT

A target-label-free, candidate-conditioned pretrained sequence-context score contains transferable information about higher-order residuals.

### H_LOW_ORDER_SHORTCUT

Any apparent sequence-model association is explained by low-order O1/O2/interaction-shape statistics, mutation order, or intrinsic single-mutation plausibility.

### H_NULL

After low-order controls, the candidate-conditioned sequence-context score has no material association with higher-order residuals.

No biological root cause may be declared from this experiment alone.

## Revealed development datasets only

1. TrpB FLIP2 two-to-many
2. GB1 Wu-2016 four-site landscape

NucB must not be downloaded, read, listed, or used.

## Target-side information contract

Allowed before evaluation:

- target WT/reference sequence;
- target candidate sequence;
- target mutation identities/positions;
- target WT/single/double measured fitness;
- frozen pretrained ESM-2 weights;
- exact low-order O1/O2/e1/e2 quantities derived only from target WT/single/double labels.

Forbidden for scoring/model inference:

- target triple fitness;
- target quad fitness;
- any target higher-order label;
- target-label fit of ESM;
- target-label calibration;
- target-label feature selection;
- target-label threshold selection;
- target-label blend weights.

Higher-order target labels are used only after all ESM scores are frozen, for evaluation.

## Immutable public sequence/model references

### GB1 full sequence

Use the 56-aa chain-A sequence from PDB 3GB1:

MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE

Four benchmark positions are 39, 40, 41, and 54.
The expected WT genotype at those positions is VDGV.

The run must abort if this invariant fails.

Public reference:
https://www.rcsb.org/structure/3GB1

### Frozen sequence model

Repository:
facebook/esm2_t6_8M_UR50D

Revision:
c731040fcd8d73dceaa04b0a8e6329b345b0f5df

Expected model.safetensors SHA256:
24c5fa474c48f3b754b86efe752d5f189d2bcd88190fa2270fc92b2ef3034189

Architecture:
EsmForMaskedLM

The model is evaluation-only, no gradients, no fine-tuning.

Frozen runtime:
- Python 3.12
- torch 2.14.0 CPU
- transformers 5.17.0
- random seed 161

## Candidate eligibility

For each dataset independently:

1. fit robust_asinh parameters using LOW-ORDER labels only;
2. derive exact transformed WT, e1, and e2;
3. a higher-order candidate is eligible only if every singleton and every pair required for exact O2 exists;
4. analyze mutation orders 3 and 4 separately.

## Label-independent sampling

For each of four cells:

- TrpB order 3
- TrpB order 4
- GB1 order 3
- GB1 order 4

select exactly 192 eligible candidates.

Selection key:

sha256("NABU_SEQ_CONTEXT_V1|161|<dataset>|<order>|<candidate_id>")

Select the 192 lexicographically smallest hashes.

Fitness labels MUST NOT participate in candidate selection.

If any cell has fewer than 192 eligible candidates, abort instead of changing sample size.

The exact selected IDs and hashes must be stored before ESM evaluation.

## Higher-order evaluation target

On the low-order-fit robust_asinh scale:

residual_HO = transformed_true_fitness - exact_O2

This is the evaluation target only.

## Frozen ESM scoring rule

For each mutation j in candidate S:

1. construct the full candidate sequence;
2. mask mutation position j while leaving all OTHER candidate mutations present;
3. compute:

candidate_preference_j =
log P(mutant_aa_j | candidate context, j masked)
-
log P(wildtype_aa_j | candidate context, j masked)

4. construct the corresponding SINGLE-mutant sequence containing only mutation j;
5. mask position j;
6. compute:

single_preference_j =
log P(mutant_aa_j | single-mutant context, j masked)
-
log P(wildtype_aa_j | single-mutant context, j masked)

7. define context shift:

context_shift_j =
candidate_preference_j - single_preference_j

Primary candidate sequence-context score:

ESM_CONTEXT_MEAN =
mean_j(context_shift_j)

Secondary stored scores:
- ESM_CONTEXT_SUM
- ESM_CONTEXT_STD
- ESM_SINGLE_MEAN
- every per-mutation candidate_preference_j
- every per-mutation single_preference_j
- every per-mutation context_shift_j

No score is selected post hoc.

## Low-order nuisance features

Use the exact 15-feature low-order vector already frozen in prior falsification work:

1. O1 / sigma_low
2. O2 / sigma_low
3. sum(e1) / sigma_low
4. mean(abs(e1)) / sigma_low
5. max(abs(e1)) / sigma_low
6. std(e1) / sigma_low
7. min(e1) / sigma_low
8. max(e1) / sigma_low
9. sum(e2) / sigma_low
10. mean(abs(e2)) / sigma_low
11. max(abs(e2)) / sigma_low
12. std(e2) / sigma_low
13. min(e2) / sigma_low
14. max(e2) / sigma_low
15. fraction(e2 > 0)

sigma_low is std of transformed LOW-ORDER labels only.

## Evaluation

Evaluation occurs only after candidate IDs and all ESM scores are frozen.

For every dataset/order cell report:

### Raw association
- Spearman(ESM_CONTEXT_MEAN, residual_HO)
- Pearson
- RMSE after source-free best affine fit, marked evaluation-only and NOT used as a gate

### Intrinsic single-mutation control
- Spearman(ESM_SINGLE_MEAN, residual_HO)

### Incremental low-order-controlled association

Evaluation-only OLS residualization:

- regress residual_HO on intercept + frozen 15 low-order nuisance features;
- regress ESM_CONTEXT_MEAN on intercept + the same nuisance features;
- regress ESM_SINGLE_MEAN on intercept + the same nuisance features.

Report:

PARTIAL_CONTEXT_SPEARMAN =
Spearman(residualized residual_HO, residualized ESM_CONTEXT_MEAN)

PARTIAL_SINGLE_SPEARMAN =
Spearman(residualized residual_HO, residualized ESM_SINGLE_MEAN)

This residualization is evaluation-only and cannot feed a future model.

### Permutation null

Within each dataset/order cell:

- permute ESM_CONTEXT_MEAN across the 192 selected candidates;
- 1000 permutations;
- seed 161;
- compute partial Spearman after the same nuisance residualization.

Report:
- empirical two-sided p-value;
- null mean/std;
- null 95th and 99th percentiles of absolute partial Spearman.

## Primary gate

The experiment declares:

PRETRAINED_SEQUENCE_CONTEXT_SIGNAL_SUPPORTED_FOR_QUADS

only if ALL conditions below hold for BOTH TrpB order-4 and GB1 order-4:

1. abs(PARTIAL_CONTEXT_SPEARMAN) >= 0.15
2. empirical permutation p <= 0.01
3. abs(PARTIAL_CONTEXT_SPEARMAN)
   >= abs(PARTIAL_SINGLE_SPEARMAN) + 0.05
4. the sign of PARTIAL_CONTEXT_SPEARMAN is the same in TrpB and GB1
5. all model/source/hash/determinism checks pass

Triples are reported as secondary consistency evidence and do not rescue a failed quad gate.

## Rejection / unresolved rules

Declare:

PRETRAINED_SEQUENCE_CONTEXT_SCALAR_REJECTED

only if BOTH quad cells satisfy:

- abs(PARTIAL_CONTEXT_SPEARMAN) < 0.08
  OR empirical p > 0.10

Otherwise, if the support gate does not pass:

PRETRAINED_SEQUENCE_CONTEXT_SCALAR_UNRESOLVED

A rejection applies ONLY to this frozen scalar candidate-conditioned ESM-2 scoring rule.
It does not reject every possible sequence-aware representation.

## Stopping rule

If SUPPORTED:
- stop diagnostics for this question;
- do not consume NucB;
- preregister exactly one architecture that uses the supported sequence-context signal;
- no architecture may be implemented before its gate is frozen.

If REJECTED:
- do not build an architecture from this scalar ESM score;
- do not tune ESM size, score weights, or thresholds on these same results;
- any richer sequence-representation test must be separately justified as a different hypothesis class.

If UNRESOLVED:
- do not build an architecture;
- preserve all results;
- root-cause status remains ROOT_CAUSE_NOT_YET_IDENTIFIED.

## Determinism

Run ESM scoring twice on the first 16 selected candidates in each cell.

Maximum absolute difference in every stored probability-derived score must be exactly 0.0 within the same runtime.

## Required outputs

- SELECTED_CANDIDATES.csv
- ESM_MUTATION_SCORES.csv
- ESM_CANDIDATE_SCORES.csv
- CELL_METRICS.json
- DECISION.json
- LOW_ORDER_DIAGNOSTICS.json
- MODEL_MANIFEST.json
- SOURCE_MANIFEST.json
- RUNTIME_MANIFEST.json
- OUTPUT_HASHES.json
- raw permutation summaries for every cell

All failed and intermediate evidence is retained.
