# NABU Frozen-Embedding Transfer Falsification 03

Date: 2026-10-01
Branch: exp/nabu-frozen-embedding-transfer-falsification
Status: preregistered discriminative test
NucB: FORBIDDEN
Phase 3: CLOSED
Governance: SCIENTIFIC_CONSTRAINTS.md on main

## Scientific question

Does a frozen pretrained protein representation contain transferable beyond-pair context information that can be decoded from one revealed protein landscape and predict higher-order residuals in another protein landscape without using any target higher-order labels for fitting?

This tests a genuinely richer representation class than Experiments 01 and 02.

It is NOT an architecture-development experiment.

## Why this is not a repair of rejected likelihood scores

Experiments 01 and 02 tested fixed scalar functions of ESM masked-likelihood outputs.

Both were rejected.

This experiment:
- does not use either rejected scalar as the representation;
- does not tune likelihood weights;
- does not change ESM model size;
- uses the frozen last hidden representation vector;
- constructs an exact vector beyond-pair residual by the same O2 algebra;
- tests cross-landscape transfer with a fixed source-only decoder.

## Competing hypotheses

### H_TRANSFERABLE_REPRESENTATION

The frozen ESM hidden representation contains beyond-pair information that transfers across proteins and adds predictive value beyond low-order biological features.

### H_LOW_ORDER_ONLY

Any cross-landscape prediction is explained by the 15 frozen low-order biological features; the pretrained representation adds no material value.

### H_SOURCE_OVERFIT_NULL

A source-trained decoder does not transfer; any apparent source fit is landscape-specific or compatible with shuffled-label control.

No biological root cause may be declared from this experiment alone.

## Development landscapes

1. TrpB FLIP2 two-to-many
2. GB1 Wu-2016 four-site landscape

NucB is forbidden.

## Fresh label-independent holdout

Use the same eligibility and hash ordering as Experiments 01 and 02:

sha256("NABU_SEQ_CONTEXT_V1|161|<dataset>|<order>|<candidate_id>")

Prior use:
- Experiment 01: ranks 0..191
- Experiment 02: ranks 192..383

This experiment MUST use:
- ranks 384..575 inclusive

for each of:
- TrpB order 3
- TrpB order 4
- GB1 order 3
- GB1 order 4

Exactly 192 candidates per cell.

No fitness label may participate in selection.

## Information contract

### Target landscape allowed before target reveal

- target WT/reference sequence;
- target candidate sequence;
- mutation identities/positions;
- target WT/single/double measured fitness;
- frozen ESM weights;
- target low-order O1/O2/e1/e2 and frozen 15-feature vector;
- source-landscape higher-order labels, only when that landscape is the declared source direction.

### Forbidden before target prediction freeze

- target triple fitness;
- target quad fitness;
- target higher-order residual;
- target-label feature selection;
- target-label standardization;
- target-label decoder fitting;
- target-label hyperparameter tuning;
- target-label early stopping;
- target-label threshold or blend selection.

## Cross-landscape directions

Two independent sealed directions:

1. GB1 source -> TrpB target
2. TrpB source -> GB1 target

For each direction:
- source higher-order labels are allowed for decoder fitting;
- target higher-order labels remain sealed until all target predictions are frozen.

## Frozen ESM model

Model:
facebook/esm2_t6_8M_UR50D

Revision:
c731040fcd8d73dceaa04b0a8e6329b345b0f5df

Expected model.safetensors:
- size: 31,384,292 bytes
- SHA256: 24c5fa474c48f3b754b86efe752d5f189d2bcd88190fa2270fc92b2ef3034189

Architecture:
EsmForMaskedLM, used in eval mode with output_hidden_states=True.

Representation:
last hidden layer, dimension expected from pinned model config.

No fine-tuning.
No gradients.
No model-size search.

## Candidate-side representation

For candidate mutation set S and any subset T among:
- empty set,
- every singleton in S,
- every pair in S,
- full S

construct sequence x_T.

Let P(S) be the full candidate mutation positions.

Define:

E_S(T) =
mean last-layer hidden vector over positions in P(S)
when the input sequence is x_T.

All E_S(T) vectors for a candidate use the SAME position set P(S).

Define vector low-order terms:

E0 = E_S(empty)

E1(j) = E_S({j}) - E0

E2(j,k) =
E_S({j,k}) - E0 - E1(j) - E1(k)

Vector O2 reconstruction:

E_O2(S) =
E0
+ sum_j E1(j)
+ sum_{j<k} E2(j,k)

Primary vector representation:

R_HO(S) =
E_S(S) - E_O2(S)

This is the frozen pretrained representation residual beyond pair order.

No target labels are used to construct R_HO.

## Biological target

Fit robust_asinh using target LOW-ORDER labels only.

Derive exact transformed target WT/e1/e2.

For each eligible higher-order candidate:

BIO_HO =
transformed_true_fitness - exact_biological_O2

BIO_HO is visible only for the declared source landscape during decoder fitting and only after prediction freeze for the target landscape.

## Low-order feature arm

Use the same frozen 15-feature low-order vector from Experiments 01/02.

## Fixed decoder

For each direction train four fixed arms on the SOURCE higher-order rows only.

Source training rows:
- 192 order-3
- 192 order-4
- total 384

No order indicator is supplied.

### LOW_ONLY

Features:
15 frozen low-order biological features.

Pipeline:
StandardScaler fit on SOURCE only
-> Ridge(alpha=1.0, solver="svd", fit_intercept=True)

### REPRESENTATION_ONLY

Features:
R_HO vector only.

Pipeline:
StandardScaler fit on SOURCE only
-> Ridge(alpha=1.0, solver="svd", fit_intercept=True)

### FULL

Features:
15 low-order features concatenated with R_HO.

Pipeline:
StandardScaler fit on SOURCE only
-> Ridge(alpha=1.0, solver="svd", fit_intercept=True)

### SHUFFLED_FULL

Same FULL features and pipeline.

SOURCE BIO_HO labels are permuted once with seed 161 before fitting.

No hyperparameter search.

## Prediction freeze

Each direction must run in a job that has:
- frozen representations for both landscapes;
- source labels for only the declared source landscape;
- NO target higher-order labels.

Freeze predictions for every target candidate before target labels become available.

Store:
- all four arm predictions;
- scaler parameters;
- Ridge coefficients/intercept;
- source row IDs;
- target row IDs;
- feature dimension;
- prediction hashes.

## Evaluation

After both direction predictions are frozen, target labels may be revealed.

For each direction report on target:

### BIO_HO residual prediction
Whole target:
- Spearman
- Pearson
- RMSE

Order 3:
- Spearman
- Pearson
- RMSE

Order 4:
- Spearman
- Pearson
- RMSE

for:
- LOW_ONLY
- REPRESENTATION_ONLY
- FULL
- SHUFFLED_FULL

### Reconstructed transformed fitness

For each arm:

predicted_transformed_fitness =
exact_target_O2 + predicted_BIO_HO

Report whole/order-3/order-4:
- Spearman versus true transformed fitness

Also report exact O2 alone.

## Primary gate

Declare:

FROZEN_EMBEDDING_TRANSFER_SUPPORTED

only if BOTH directions satisfy ALL on target order-4:

1. FULL BIO_HO Spearman >= 0.25
2. FULL BIO_HO Spearman >= LOW_ONLY BIO_HO Spearman + 0.10
3. FULL BIO_HO Spearman >= SHUFFLED_FULL BIO_HO Spearman + 0.15
4. REPRESENTATION_ONLY BIO_HO Spearman >= 0.15
5. FULL reconstructed-fitness Spearman >= exact-O2 reconstructed-fitness Spearman + 0.10
6. order-3 FULL BIO_HO Spearman >= order-3 LOW_ONLY BIO_HO Spearman - 0.02
7. all leakage/source/model/hash/determinism checks pass

Both directions must pass.

## Rejection / unresolved

Declare:

FROZEN_EMBEDDING_TRANSFER_REJECTED

only if BOTH directions satisfy:

- order-4 FULL BIO_HO Spearman < 0.10
AND
- order-4 FULL minus LOW_ONLY BIO_HO Spearman < 0.05

Otherwise, if the support gate does not pass:

FROZEN_EMBEDDING_TRANSFER_UNRESOLVED

## Stopping rule

If SUPPORTED:
- stop representation diagnostics;
- do not consume NucB;
- preregister exactly one architecture using the supported representation prior.

If REJECTED:
- do not build an architecture from this ESM-2 8M hidden-representation transfer path;
- do not tune Ridge alpha, order weights, feature pooling, or layer choice on these revealed results;
- do not test another ESM-2 layer as a repair;
- any next prior class must be genuinely different or require additional measurements.

If UNRESOLVED:
- do not build an architecture;
- preserve all outputs;
- define one new discriminative question before further work.

## Determinism

Representation extraction:
- replay the identical full sequence/subset stream;
- maximum absolute representation difference must be exactly 0.0.

Decoder:
- fit each source arm twice with identical inputs;
- maximum absolute target prediction difference must be exactly 0.0.

## Required outputs

Preparation:
- SELECTED_CANDIDATES.csv
- LOW_ORDER_DIAGNOSTICS.json
- SOURCE_MANIFEST.json
- source-label artifacts split by landscape
- target labels never co-resident with a direction prediction job

Representation:
- SUBSET_REPRESENTATIONS.csv or equivalent lossless numeric artifact
- CANDIDATE_R_HO.csv
- MODEL_MANIFEST.json
- RUNTIME_MANIFEST.json
- REPRESENTATION_DIAGNOSTICS.json
- hashes

Direction prediction bundles:
- GB1_TO_TRPB_PREDICTIONS.csv
- TRPB_TO_GB1_PREDICTIONS.csv
- decoder coefficients/scalers/manifests
- source/target ID hashes
- prediction hashes
- determinism checks

Final evaluation:
- TRANSFER_METRICS.json
- DECISION.json
- JOINED_EVIDENCE.csv
- OUTPUT_HASHES.json

All intermediate, negative, and failed evidence is retained.
