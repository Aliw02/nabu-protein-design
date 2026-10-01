# Evidence Record — Sequence-Context Falsification 01

Canonical experiment branch:
`exp/nabu-sequence-context-falsification`

Canonical result commit:
`c668ab156dbc22d02e29e5648fd0765938497d31`

Canonical scientific run:
`36875022535`

Complete evidence artifact:
- ID: `11168668376`
- digest: `sha256:66bd9258f4879491d61cc45713d6b54a1f609a874fcc85e19bd64b793b67d9b4`

## Decision

`PRETRAINED_SEQUENCE_CONTEXT_SCALAR_REJECTED`

`ROOT_CAUSE_NOT_YET_IDENTIFIED`

NucB consumed: false  
Phase 3 opened: false

## Primary quad result

TrpB order 4:
- n = 192
- raw context Spearman = -0.10087011366410764
- raw context Pearson = -0.10885684776338019
- partial context Spearman after frozen low-order nuisance control = **-0.028251566611507472**
- partial single-mutation control Spearman = 0.048849795187586476
- empirical permutation p = **0.7202797202797203**
- permutation null |rho| p95 = 0.15285635325394029
- permutation null |rho| p99 = 0.19092228318368012
- all three preregistered quad support conditions = FAIL

GB1 order 4:
- n = 192
- raw context Spearman = -0.22074980332582808
- raw context Pearson = -0.16530090124991553
- partial context Spearman after frozen low-order nuisance control = **-0.03728501478447224**
- partial single-mutation control Spearman = 0.009202994872907795
- empirical permutation p = **0.6113886113886113**
- permutation null |rho| p95 = 0.14539587458969697
- permutation null |rho| p99 = 0.1871617204242737
- all three preregistered quad support conditions = FAIL

Both quad cells satisfy the preregistered weak-rejection condition.

## Secondary triple result

TrpB order 3:
- n = 192
- raw context Spearman = -0.10603450614437242
- partial context Spearman = 0.07848316740362965
- partial single Spearman = -0.1309816075739902
- empirical permutation p = 0.3016983016983017

GB1 order 3:
- n = 192
- raw context Spearman = -0.23622263516262923
- partial context Spearman = -0.006227450288907577
- partial single Spearman = 0.01652571684344736
- empirical permutation p = 0.9310689310689311

## Frozen ESM identity

- model: `facebook/esm2_t6_8M_UR50D`
- revision: `c731040fcd8d73dceaa04b0a8e6329b345b0f5df`
- model.safetensors SHA256: `24c5fa474c48f3b754b86efe752d5f189d2bcd88190fa2270fc92b2ef3034189`
- architecture: EsmForMaskedLM
- fine-tuned: false
- deterministic replay max absolute difference: **0.0**

## Runtime

- Python 3.12.14
- numpy 2.4.6
- pandas 3.0.6
- torch 2.14.0+cpu
- transformers 5.17.0
- huggingface_hub 1.33.0
- safetensors 0.8.0
- seed 161

## Sampling and scoring volume

Four label-independent cells:
- TrpB triples: 192
- TrpB quads: 192
- GB1 triples: 192
- GB1 quads: 192
- total selected candidates: 768

Scoring:
- candidate masked requests: 2,688
- unique single-mutant masked requests: 295
- total masked forward examples: 2,983
- mutation-level score rows: 2,688
- candidate-level score rows: 768
- unique TrpB mutations represented: 220
- unique GB1 mutations represented: 75

## Low-order support

TrpB:
- fit = 8,633
- validation = 2,158
- test = 217,507
- fit WT/singles/doubles = 1 / 316 / 8,316
- test triples/quads = 87,989 / 129,518
- exact-O2 eligible triples/quads = 24,695 / 14,816
- e1/e2 count = 316 / 5,670
- sigma_low = 1.3851019211772997
- robust_asinh center = -0.0007551379217456
- robust_asinh scale = 0.0368309720526426

GB1:
- fit = 2,168
- test = 147,193
- fit WT/singles/doubles = 1 / 76 / 2,091
- test triples/quads = 26,019 / 121,174
- exact-O2 eligible triples/quads = 25,035 / 112,286
- e1/e2 count = 76 / 2,091
- sigma_low = 2.187748986028811
- robust_asinh center = 0.0196384306666
- robust_asinh scale = 0.0196384306666

## Source identity

TrpB:
- MD5: `a611408d2db06907a023ddc8a1d94c12`
- SHA256: `2676d06452f4f5dbd52cb5c8cad17d34c94ac8307bda581b527be208a3925ba5`

GB1:
- size = 2,993,024 bytes
- Git blob SHA = `a1e9f5146a9441210eaf88accfa660f02d80cbb7`
- SHA256 = `02ef9bc1784a67ee14cf4024596a8fcfc7129779911c6589ed896373552125b1`
- PDB reference = 3GB1
- positions = 39, 40, 41, 54
- WT genotype = VDGV

## Artifact identities

Prediction-visible:
- artifact ID: 11169055105
- digest: `sha256:2b823698b53dec82e4a05f0c5de99b8e7203e005253f339f9564ce25c60b8b34`

Frozen scores:
- artifact ID: 11168806614
- digest: `sha256:3abcfe0b0d279f038096d330c0dd9ad998ec64af61060326efaa0f3b1b89b25e`

Sealed reveal:
- artifact ID: 11168731241
- digest: `sha256:b52d9091fc05d64439cb22f36fd7b08f9acae9c6f261ef2c178c3b46498da3df`

Complete evidence:
- artifact ID: 11168668376
- digest: `sha256:66bd9258f4879491d61cc45713d6b54a1f609a874fcc85e19bd64b793b67d9b4`

## Canonical file hashes

- SELECTED_CANDIDATES.csv: `2749db7b5f25ded4007b8c26cce59cf09083eff888af797294ec7f159f8ebfaf`
- LOW_ORDER_DIAGNOSTICS.json: `fe425036e9d6381430e5b089b8340f80e74a200f467170b3a345a494039f5f86`
- SOURCE_MANIFEST.json: `8a99d8b0bd620e38943eefab6fdfd1473b24e522c2e187ba9e02bade018b5d08`
- ESM_CANDIDATE_SCORES.csv: `18a41f56c6f835e4f4f6fce44dac85fe8a4532120ebbc9d3815b96e2530e3a3b`
- ESM_MUTATION_SCORES.csv: `46d1a082044188ceea09d9b45458d13007c885586e312c1b7e9c763810c6a624`
- MODEL_MANIFEST.json: `08ecfecdb3b2ced90394980f3877306dd4e2e0c2a3a9020e4597e38873d12d99`
- RUNTIME_MANIFEST.json: `1771995bef5bdbc7bcb7fdd13e76b324f82aa111e4ee418fd1b2969852173390`
- SCORING_DIAGNOSTICS.json: `11506d3459c6c87d8f1dc382f0e096013dc7df27104e8be54f28f511bc542826`
- REVEAL_TARGETS.csv: `43a97dcda66bc9ec8a708dd4babcd3b290103137e82a5c527006193e516ae26d`
- CELL_METRICS.json: `42a6dbe56e9aab97d511057ce6bab41a465acddd4ae235dbcd3696e07c2f5960`
- DECISION.json: `c09b18a794b8af1140004fdf171d2cec73fc019ef1f34a6825d230848adcdced`
- EVALUATION_DIAGNOSTICS.json: `12521cc1bd2001fb95a10faf3e6dcb666a88bad8c6a22044a3964ead15f2150b`
- JOINED_EVIDENCE.csv: `499fc6422b89114223365273dca83c7b0772f123ff4b4adae4cf8053b06a122d`
- PERMUTATION_NULL.csv: `859f24a8eb003fe2c6e5c8d94bdb4b3c2f87e1dc5574902c1d9da50e7bcfe4a2`

## Interpretation

The raw candidate-conditioned ESM context score had modest correlations in some cells, especially GB1, but those associations disappeared after controlling for the frozen low-order feature set.

The primary quad partial correlations were near zero and well inside the permutation null.

Therefore this scalar candidate-conditioned ESM preference-shift is not supported as the missing transferable higher-order signal and must not be turned into an architecture.

A richer frozen sequence representation remains a distinct, unresolved hypothesis class. It cannot be treated as a patch of this score and requires a new preregistered discriminative protocol.
