# NABU Pretrained Möbius Falsification 02

Date: 2026-10-01
Branch: exp/nabu-pretrained-mobius-falsification
Status: preregistered discriminative test
NucB: FORBIDDEN
Phase 3: CLOSED
Governance: SCIENTIFIC_CONSTRAINTS.md on main

## Scientific question

Does a frozen pretrained protein language model contain explicit beyond-pair non-additivity that tracks the biological higher-order residual y(S)-O2, without using any target higher-order labels for scoring or fitting?

This is a falsification experiment, NOT an architecture-development experiment.

## Why this is a distinct hypothesis class

Sequence-Context Falsification 01 tested one scalar:
the mean change in masked mutant-vs-WT preference between candidate context and single-mutant context.

That scalar was rejected.

This experiment does NOT tune that score.

Instead it constructs a model-side higher-order residual using the same algebraic decomposition as the biological target:

MODEL_HO(S) = MODEL_F(S) - MODEL_O2(S)

where MODEL_O2 is reconstructed exactly from model-side singleton and pair terms.

Thus the tested hypothesis is explicit pretrained-model non-additivity beyond pair order.

## Competing hypotheses

### H_PRETRAINED_MOBIUS

Frozen pretrained model beyond-pair non-additivity contains transferable information about biological higher-order residuals after low-order nuisance control.

### H_LOWER_ORDER_CONTEXT

Any apparent association is already explained by low-order biology or by the simpler scalar context-shift score rejected in Experiment 01.

### H_NULL

Frozen pretrained model beyond-pair non-additivity has no material association with biological higher-order residuals after controls.

No biological root cause may be declared from this experiment alone.

## Revealed development datasets

1. TrpB FLIP2 two-to-many
2. GB1 Wu-2016 four-site landscape

NucB must not be downloaded, read, listed, or used.

## Fresh label-independent sample

Reuse the exact immutable source identities and low-order eligibility rules from Sequence-Context Falsification 01.

For each of four cells:
- TrpB order 3
- TrpB order 4
- GB1 order 3
- GB1 order 4

sort eligible candidates by:

sha256("NABU_SEQ_CONTEXT_V1|161|<dataset>|<order>|<candidate_id>")

Experiment 01 used ranks 0..191.

This experiment MUST use ranks 192..383 inclusive.

Exactly 192 candidates per cell.

Fitness labels MUST NOT participate in candidate selection.

If any cell has fewer than 384 eligible candidates, abort.

Store exact IDs, hashes, and rank before any ESM scoring.

## Target-side information contract

Allowed before evaluation:
- WT/reference sequence;
- candidate sequence;
- mutation identities and positions;
- WT/single/double measured fitness;
- frozen pretrained ESM-2 weights;
- exact target low-order O1/O2/e1/e2 quantities from WT/singles/doubles only.

Forbidden before score freeze:
- target triple fitness;
- target quad fitness;
- target higher-order residuals;
- target-label model fitting;
- target-label feature selection;
- target-label calibration;
- target-label threshold selection;
- target-label blend weights.

Higher-order labels are reveal-only evaluation data.

## Frozen model

Model:
facebook/esm2_t6_8M_UR50D

Revision:
c731040fcd8d73dceaa04b0a8e6329b345b0f5df

Expected model.safetensors:
- size: 31,384,292 bytes
- SHA256: 24c5fa474c48f3b754b86efe752d5f189d2bcd88190fa2270fc92b2ef3034189

Architecture:
EsmForMaskedLM

Execution:
- eval mode
- no gradients
- no fine-tuning
- seed 161
- CPU
- Python 3.12
- torch 2.14.0
- transformers 5.17.0

No model-size search is allowed.

## Model-side scalar function

For a candidate mutation set S and any non-empty subset T subseteq S:

Construct sequence x_T containing exactly mutations in T.

For each mutation j in T:
1. mask position j in x_T while leaving all other mutations in T present;
2. compute:
   log_odds_j(T) =
   log P(mutant_aa_j | x_T with j masked)
   -
   log P(wildtype_aa_j | x_T with j masked)

Define:

MODEL_F(T) = sum over j in T of log_odds_j(T)

and MODEL_F(empty) = 0.

For candidate S:

MODEL_E1(j) = MODEL_F({j})

MODEL_E2(j,k) =
MODEL_F({j,k})
- MODEL_F({j})
- MODEL_F({k})

MODEL_O2(S) =
sum_j MODEL_E1(j)
+
sum_{j<k} MODEL_E2(j,k)

Primary pretrained beyond-pair score:

MODEL_HO(S) =
MODEL_F(S) - MODEL_O2(S)

For triples this is model-side third-order residual.

For quads this is total model-side beyond-pair residual, directly analogous to biological transformed target minus exact biological O2.

## Scalar controls

On the SAME fresh candidates also compute:

### CONTEXT_MEAN_CONTROL

The exact scalar from Experiment 01:
mean over mutations of candidate-context preference minus single-mutant-context preference.

### SINGLE_MEAN_CONTROL

Mean singleton masked mutant-vs-WT preference.

No post-hoc score selection.

## Biological evaluation target

Fit robust_asinh using LOW-ORDER target labels only.

Derive exact transformed WT/e1/e2 from target WT/singles/doubles.

For eligible candidate S:

BIO_HO(S) =
transformed_true_fitness(S)
-
exact_biological_O2(S)

Target higher-order labels are used only after all model-side scores are frozen.

## Low-order nuisance vector

Use the same frozen 15-feature low-order vector from Experiment 01:

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

sigma_low is derived from transformed LOW-ORDER labels only.

## Evaluation

For every dataset/order cell report:

### Raw
- Spearman(MODEL_HO, BIO_HO)
- Pearson(MODEL_HO, BIO_HO)
- Spearman(CONTEXT_MEAN_CONTROL, BIO_HO)
- Spearman(SINGLE_MEAN_CONTROL, BIO_HO)

### Low-order-controlled
Evaluation-only OLS residualization:

- residualize BIO_HO on intercept + frozen 15 low-order features;
- residualize MODEL_HO on the same design;
- residualize CONTEXT_MEAN_CONTROL on the same design;
- residualize SINGLE_MEAN_CONTROL on the same design.

Report:
- PARTIAL_MOBIUS_SPEARMAN
- PARTIAL_CONTEXT_SPEARMAN
- PARTIAL_SINGLE_SPEARMAN

No residualized quantity may feed a future architecture from this run.

### Permutation null

Within each cell:
- permute MODEL_HO across the 192 candidates;
- 1000 permutations;
- seed 161;
- recompute partial Spearman after the same nuisance residualization.

Report:
- empirical two-sided p;
- null mean/std;
- null absolute-rho p95;
- null absolute-rho p99.

## Primary gate

Declare:

PRETRAINED_MOBIUS_SIGNAL_SUPPORTED_FOR_QUADS

only if BOTH TrpB order-4 and GB1 order-4 satisfy ALL:

1. abs(PARTIAL_MOBIUS_SPEARMAN) >= 0.15
2. empirical permutation p <= 0.01
3. abs(PARTIAL_MOBIUS_SPEARMAN)
   >= abs(PARTIAL_CONTEXT_SPEARMAN) + 0.05
4. abs(PARTIAL_MOBIUS_SPEARMAN)
   >= abs(PARTIAL_SINGLE_SPEARMAN) + 0.05
5. PARTIAL_MOBIUS_SPEARMAN has the same sign in TrpB and GB1
6. all source/model/hash/determinism/leakage checks pass

Triples are secondary consistency evidence only and cannot rescue a failed quad gate.

## Rejection / unresolved rules

Declare:

PRETRAINED_MOBIUS_SCALAR_REJECTED

only if BOTH quad cells satisfy either:
- abs(PARTIAL_MOBIUS_SPEARMAN) < 0.08
  OR
- empirical permutation p > 0.10

Otherwise, if the support gate does not pass:

PRETRAINED_MOBIUS_SCALAR_UNRESOLVED

This rejection applies only to this exact pretrained likelihood-based Möbius scalar.

## Stopping rule

If SUPPORTED:
- stop diagnostics for this hypothesis;
- do not consume NucB;
- preregister one architecture using MODEL_HO before implementation.

If REJECTED:
- do not build an architecture from ESM likelihood algebra;
- do not tune model size, weights, score combinations, or thresholds on these revealed results;
- classify frozen ESM likelihood-scalar priors as exhausted for this development path;
- any next sequence-aware test must use a genuinely richer representation class, not likelihood-score repair.

If UNRESOLVED:
- do not build an architecture;
- preserve all outputs;
- root-cause status remains ROOT_CAUSE_NOT_YET_IDENTIFIED.

## Determinism

Replay the complete scoring calculation on the first 8 candidates in each cell.

Maximum absolute difference for:
- every masked log odds,
- every MODEL_F subset score,
- MODEL_O2,
- MODEL_HO,
- CONTEXT_MEAN_CONTROL,
- SINGLE_MEAN_CONTROL

must be exactly 0.0 within the same runtime.

## Required outputs

Prediction-visible:
- SELECTED_CANDIDATES.csv
- LOW_ORDER_DIAGNOSTICS.json
- SOURCE_MANIFEST.json
- VISIBLE_HASHES.json

Frozen scoring:
- MASKED_LOG_ODDS.csv
- SUBSET_MODEL_F.csv
- CANDIDATE_MODEL_SCORES.csv
- MODEL_MANIFEST.json
- RUNTIME_MANIFEST.json
- SCORING_DIAGNOSTICS.json
- SCORE_HASHES.json

Sealed reveal:
- REVEAL_TARGETS.csv
- REVEAL_HASHES.json

Evaluation:
- CELL_METRICS.json
- DECISION.json
- EVALUATION_DIAGNOSTICS.json
- JOINED_EVIDENCE.csv
- PERMUTATION_NULL.csv
- OUTPUT_HASHES.json

Every intermediate and failed result is retained.
