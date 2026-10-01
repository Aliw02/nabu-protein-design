# Pretrained Möbius Falsification 02 — Canonical Result

Date: 2026-10-01
Branch: exp/nabu-pretrained-mobius-falsification

## Canonical execution

Scientific run:
- GitHub Actions run: 36878521685
- head SHA: 829c082546edf5110af4b92ec2aefb7046ac7e18
- workflow conclusion: success

Operational-only failed run before the canonical run:
- 36877363477 — determinism replay implementation issue; no scientific decision was accepted from that run.
The canonical run replayed the identical request stream and passed exact determinism.

## Frozen scientific decision

PRETRAINED_MOBIUS_SCALAR_REJECTED

Next action:
DO_NOT_BUILD_FROM_ESM_LIKELIHOOD_ALGEBRA;RICHER_REPRESENTATION_REQUIRES_NEW_HYPOTHESIS

Root-cause status:
ROOT_CAUSE_NOT_YET_IDENTIFIED

NucB consumed: false
Phase 3 opened: false

This rejection applies to the exact pretrained likelihood-based Möbius scalar defined in PROTOCOL.md. It does not reject every possible frozen sequence representation.

## Fresh holdout contract

This experiment did not reuse the candidates revealed in Sequence-Context Falsification 01.

Selection:
- same immutable hash ordering rule;
- Experiment 01 used ranks 0..191;
- this experiment used ranks 192..383;
- 192 candidates in each dataset/order cell;
- total selected candidates: 768.

Integrity:
- fresh holdout ranks exact: true
- visible hashes exact: true
- frozen score hashes exact: true
- reveal hashes exact: true
- model revision exact: true
- model weights hash exact: true
- deterministic replay exact: true
- NucB absent: true
- Phase 3 closed: true

## Primary result — TrpB order 4

n = 192
selection ranks = 192..383

Raw:
- MODEL_HO Spearman: 0.000849428966714592
- MODEL_HO Pearson: 0.006316303563980276
- context-mean control Spearman: -0.07265245639258879
- context-mean control Pearson: -0.0887263356033848
- single-mean control Spearman: 0.36412527466565403
- single-mean control Pearson: 0.4187497677538507
- BIO_HO std: 2.9366571713616723
- MODEL_HO std: 0.11973782673342327

After frozen low-order nuisance residualization:
- partial Möbius Spearman: 0.07588401649350297
- partial context-control Spearman: 0.028480454656430572
- partial single-control Spearman: 0.17305462116485368
- |Möbius| - |context|: 0.0474035618370724
- |Möbius| - |single|: -0.0971706046713507

Permutation null:
- permutations: 1000
- empirical two-sided p: 0.3096903096903097
- null mean: 0.0013238138512871984
- null SD: 0.07546626672629635
- null |rho| p95: 0.14727563885196537
- null |rho| p99: 0.18760428817513491

Primary gate:
- abs(partial Möbius) >= 0.15: FAIL
- permutation p <= 0.01: FAIL
- beats context by >= 0.05: FAIL
- beats single by >= 0.05: FAIL

Weak-rejection condition: true

## Primary result — GB1 order 4

n = 192
selection ranks = 192..383

Raw:
- MODEL_HO Spearman: -0.2251054580473646
- MODEL_HO Pearson: -0.20343840185773593
- context-mean control Spearman: -0.3764564061525107
- context-mean control Pearson: -0.37469584612717416
- single-mean control Spearman: 0.19752529636763153
- single-mean control Pearson: 0.1963974197375997
- BIO_HO std: 3.0267788808807703
- MODEL_HO std: 0.5909998939466803

After frozen low-order nuisance residualization:
- partial Möbius Spearman: -0.0046879662534248445
- partial context-control Spearman: -0.03593372758592627
- partial single-control Spearman: -0.07710305726609339
- |Möbius| - |context|: -0.031245761332501426
- |Möbius| - |single|: -0.07241509101266855

Permutation null:
- permutations: 1000
- empirical two-sided p: 0.9590409590409591
- null mean: -0.00048108537015435516
- null SD: 0.07296161556479079
- null |rho| p95: 0.14242567072674497
- null |rho| p99: 0.19391200187179555

Primary gate:
- abs(partial Möbius) >= 0.15: FAIL
- permutation p <= 0.01: FAIL
- beats context by >= 0.05: FAIL
- beats single by >= 0.05: FAIL

Weak-rejection condition: true

Cross-landscape quad sign consistency:
- TrpB partial Möbius: +0.07588401649350297
- GB1 partial Möbius: -0.0046879662534248445
- same sign: false

## Secondary result — TrpB order 3

n = 192

Raw:
- Möbius Spearman: -0.1008548544611128
- Möbius Pearson: -0.039523834307967735
- context control Spearman: -0.047932547540894675
- single control Spearman: 0.19685558690285654
- BIO_HO std: 2.055701649751193
- MODEL_HO std: 0.438504255313117

Controlled:
- partial Möbius Spearman: 0.11466273770447334
- partial context Spearman: 0.11814692238830264
- partial single Spearman: -0.20688936060548518

Permutation:
- p: 0.12087912087912088
- null |rho| p95: 0.14323517144562298
- null |rho| p99: 0.17638819751512355

## Secondary result — GB1 order 3

n = 192

Raw:
- Möbius Spearman: -0.05345976995903753
- Möbius Pearson: -0.20705283743639324
- context control Spearman: -0.1761081572308277
- single control Spearman: 0.0006900550687681416
- BIO_HO std: 1.6035750371493689
- MODEL_HO std: 0.33442231288898566

Controlled:
- partial Möbius Spearman: -0.11310799446599573
- partial context Spearman: 0.052878224778232924
- partial single Spearman: -0.0759297941024876

Permutation:
- p: 0.13686313686313686
- null |rho| p95: 0.14476702587960827
- null |rho| p99: 0.18301947413395545

## Frozen model identity

Model:
facebook/esm2_t6_8M_UR50D

Revision:
c731040fcd8d73dceaa04b0a8e6329b345b0f5df

model.safetensors:
- size: 31,384,292 bytes
- SHA256: 24c5fa474c48f3b754b86efe752d5f189d2bcd88190fa2270fc92b2ef3034189

Runtime:
- Python: 3.12.14
- numpy: 2.4.6
- pandas: 3.0.6
- torch: 2.14.0+cpu
- transformers: 5.17.0
- huggingface_hub: 1.33.0
- safetensors: 0.8.0
- seed: 161
- batch size: 16
- torch threads: 2

## Scoring volume

- selected candidates: 768
- unique masked requests: 7,626
- masked score rows: 7,626
- subset MODEL_F rows: 6,912
- candidate score rows: 768
- GB1 masked requests: 3,771
- TrpB masked requests: 3,855

Determinism replay:
- candidates: 32
- request stream: 7,626
- max masked-score difference: 0.0
- max derived-score difference: 0.0
- overall max difference: 0.0

## Immutable artifact identities

Complete evidence:
- artifact ID: 11171700565
- digest: sha256:a8cbdd632cd3acff8e18ef72e351dfc214cda1ea2bec589325120f90abb5b110

Frozen scores:
- artifact ID: 11170986879
- digest: sha256:dd3d371e81ad66e10b2ac37f1a2600ec2aae3f31dff31d08afd6c8507dd2612b

Prediction-visible:
- artifact ID: 11170650128
- digest: sha256:2f3097776a4781537c355ec6365eb3e6d9c8159f11d52a78f01312c5a49c9006

Sealed reveal:
- artifact ID: 11169958947
- digest: sha256:0f9c6472370a9980daaa10b48af9f146c957559d5299fbb2424aefeae2b3ed40

## Canonical file hashes

Prediction-visible:
- SELECTED_CANDIDATES.csv: 093f6399deb9a5a1d556199dad31dec1a50625546ef4b3646c2cef9f1ece7869
- LOW_ORDER_DIAGNOSTICS.json: fe425036e9d6381430e5b089b8340f80e74a200f467170b3a345a494039f5f86
- SOURCE_MANIFEST.json: bba564691d778a7f36d5827a8d805ab8e3bea76a6689a6df01b25a41f4140ac1

Frozen scores:
- CANDIDATE_MODEL_SCORES.csv: 03d0f8499b0bfc927d695f8ef3c56ccd235c4a7565ed66a58749c380492ed938
- MASKED_LOG_ODDS.csv: 391e68552acec0ef8f7db40176cffda31246a23620f5ed1ed7a5679694409200
- MODEL_MANIFEST.json: bffc9e1390617449b919ec8c1d591e13b3ab0269f5bb9e3eb861fa3293975d68
- RUNTIME_MANIFEST.json: 4b052685cb32b08e68e4b677fb39cd58b79ef4417464f3aa6f9359fb15ceb775
- SCORING_DIAGNOSTICS.json: 2b6b462fd8661afb6aac097bb868a1c8663ad6fbcff945525b6b07cba689a858
- SUBSET_MODEL_F.csv: 5da266d5c877d96801e8c82fed0555d665e3c93e661f665845dd812dc97cc1f1

Sealed reveal:
- REVEAL_TARGETS.csv: 69e6f57a517bbfd68a759dc73046ad2b2dfeb8d04cd7befd41bb8520f28f2bac

Evaluation:
- CELL_METRICS.json: bfe77fafd26e2569680837c79ba230fea938bfcdb498681bbf9d83c174d8ad66
- DECISION.json: 38fca188f9875419816d5960f25e828c01fbd2e4644425bf6d43bb33230264e1
- EVALUATION_DIAGNOSTICS.json: 4403a9265afcc3c4df311afd5a60bdbe82fded11e0f0db0fae5469dababde26d
- JOINED_EVIDENCE.csv: 87914d162c6bd00401bd9156bb1e7cbdc2bf8acacbcae76affa74c7d0b2f459c
- PERMUTATION_NULL.csv: e661428a1a7fdc4c03b1808dc8441dd5da491ec96d08dca634b6534ff8add0f2

## Interpretation boundary

Experiment 01 rejected one candidate-conditioned ESM preference-shift scalar.

Experiment 02 used a fresh, disjoint holdout and tested a stronger likelihood-algebra hypothesis: explicit pretrained-model beyond-pair non-additivity constructed with the same O2/Möbius decomposition used for the biological residual.

That hypothesis is also rejected by the preregistered quad rule.

The strongest evidence is after low-order nuisance removal:
- TrpB quad partial Möbius rho = +0.075884, p = 0.30969
- GB1 quad partial Möbius rho = -0.004688, p = 0.95904

The signs disagree across landscapes and neither is outside the permutation null.

Therefore:
- do not build an architecture from this MODEL_HO scalar;
- do not tune ESM size, weights, combinations, or thresholds on these revealed results;
- retire ESM likelihood-scalar priors for this development path;
- do not consume NucB;
- Phase 3 remains closed;
- a richer frozen representation is a separate hypothesis class and requires a new preregistered discriminative protocol.
