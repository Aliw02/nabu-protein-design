# Sequence-Context Falsification 01 — Canonical Result

Date: 2026-10-01
Branch: exp/nabu-sequence-context-falsification

## Canonical execution

Scientific run:
- GitHub Actions run: 36875022535
- head SHA: c09667948bbf1b9912322175ae6c0d50fcbab98f
- workflow conclusion: success

Operational-only failed runs before the canonical run:
- 36874732209 — missing legacy requirements file on clean main-derived branch
- 36874811622 — JSON diagnostic-print serialization bug after preparation
Neither reached scientific evaluation and neither is a scientific result.

## Frozen scientific decision

PRETRAINED_SEQUENCE_CONTEXT_SCALAR_REJECTED

Next action:
DO_NOT_BUILD_FROM_THIS_SCALAR_SCORE

Root-cause status:
ROOT_CAUSE_NOT_YET_IDENTIFIED

NucB consumed: false
Phase 3 opened: false

This rejection applies only to the preregistered scalar candidate-conditioned ESM-2 context-preference-shift score. It does not reject every possible sequence-aware representation.

## Integrity checks

All passed:

- model revision exact: true
- model weights hash exact: true
- deterministic replay exact: true
- NucB not consumed: true
- Phase 3 closed: true
- selected row count exact: true
- score row count exact: true
- reveal row count exact: true

ESM deterministic replay:
- candidates replayed: 64
- mutation requests replayed: 224
- maximum absolute probability-score difference: 0.0

## Frozen model

Model:
facebook/esm2_t6_8M_UR50D

Revision:
c731040fcd8d73dceaa04b0a8e6329b345b0f5df

model.safetensors:
- size: 31,384,292 bytes
- SHA256: 24c5fa474c48f3b754b86efe752d5f189d2bcd88190fa2270fc92b2ef3034189

Other model files:
- config.json SHA256: facb9355e8252149f6f99dff5fd1e32e890c3c655362b261d0c81cd6b5839c85
- special_tokens_map.json SHA256: 3aedcd4211c0d43aec4e607ff60a63255f3174ead795e997350f09a5f8cd9ee1
- tokenizer_config.json SHA256: 7e9161ecdb548ec45a41cbc6b24aa4476fdd418461f491c4207baa99419a29ad
- vocab.txt SHA256: 0b82cc0a7c7cf9e567b1e5892d793285b9fbae822c964ca48696f7db44598e03

Model execution:
- architecture: EsmForMaskedLM
- eval mode: true
- fine-tuned: false
- gradients: disabled
- max position embeddings: 1026
- maximum sequence length observed: 389

## Runtime

- Python: 3.12.14
- Linux: 6.17.0-1022-azure x86_64
- seed: 161
- torch threads: 2
- numpy: 2.4.6
- pandas: 3.0.6
- torch: 2.14.0+cpu
- transformers: 5.17.0
- huggingface_hub: 1.33.0
- safetensors: 0.8.0

## Source identity and low-order state

### TrpB

Source:
- MD5: a611408d2db06907a023ddc8a1d94c12
- SHA256: 2676d06452f4f5dbd52cb5c8cad17d34c94ac8307bda581b527be208a3925ba5

Counts:
- fit: 8,633
- validation: 2,158
- test: 217,507
- fit WT: 1
- fit singles: 316
- fit doubles: 8,316
- test triples: 87,989
- test quads: 129,518

Exact low-order support:
- eligible triples: 24,695
- eligible quads: 14,816
- selected triples: 192
- selected quads: 192

Reference:
- length: 389
- SHA256: 8a0792def4d2f8b9fa91eb32613ade22e1db02a19372e44363631d22b33e979d

Transformed low-order state:
- e1 count: 316
- e2 count: 5,670
- sigma_low: 1.3851019211772997
- robust_asinh center: -0.0007551379217456
- robust_asinh scale: 0.0368309720526426
- scale source: MAD

### GB1

Source:
- size: 2,993,024 bytes
- Git blob SHA: a1e9f5146a9441210eaf88accfa660f02d80cbb7
- SHA256: 02ef9bc1784a67ee14cf4024596a8fcfc7129779911c6589ed896373552125b1
- rows: 149,361

Counts:
- fit: 2,168
- fit WT: 1
- fit singles: 76
- fit doubles: 2,091
- test: 147,193
- test triples: 26,019
- test quads: 121,174

Exact low-order support:
- eligible triples: 25,035
- eligible quads: 112,286
- selected triples: 192
- selected quads: 192

Reference:
- PDB: 3GB1
- length: 56
- positions: 39, 40, 41, 54
- WT genotype: VDGV
- reference SHA256: 7e859d82171047700fd3e9632f7a47eab4a39baedc8c3316d2fc62d3ce2260bb

Transformed low-order state:
- e1 count: 76
- e2 count: 2,091
- sigma_low: 2.187748986028811
- robust_asinh center: 0.0196384306666
- robust_asinh scale: 0.0196384306666
- scale source: MAD

## ESM scoring volume

- selected candidates: 768
- candidate masked requests: 2,688
- unique single-mutant masked requests: 295
- total masked forward examples: 2,983
- mutation-level score rows: 2,688
- candidate-level score rows: 768
- unique TrpB mutations represented: 220
- unique GB1 mutations represented: 75
- TrpB sequence length: 389
- GB1 sequence length: 56

## Results — TrpB order 3

n = 192

Raw candidate-conditioned context score:
- Spearman: -0.10603450614437242
- Pearson: -0.0904677444589181

Evaluation-only affine fit:
- intercept: -1.6317080352094178
- slope: -0.42480081930529784
- RMSE: 2.2397792115516264
- target SD: 2.2490014983816726
- normalized RMSE: 0.9958993860890345

Intrinsic single-mutation control:
- Spearman: 0.21012939804139658
- Pearson: 0.21158686734704785

After low-order nuisance residualization:
- partial context Spearman: 0.07848316740362965
- partial single Spearman: -0.1309816075739902
- |context| - |single|: -0.05249844017036054

Permutation null:
- permutations: 1,000
- empirical two-sided p: 0.3016983016983017
- null mean: 0.0003145972926782952
- null SD: 0.07403353193375374
- null |rho| p95: 0.14258317961099204
- null |rho| p99: 0.18904548598866078

## Results — TrpB order 4

n = 192

Raw candidate-conditioned context score:
- Spearman: -0.10087011366410764
- Pearson: -0.10885684776338019

Evaluation-only affine fit:
- intercept: -2.575326157538204
- slope: -0.7964171869070158
- RMSE: 2.749188048987481
- target SD: 2.765622939431146
- normalized RMSE: 0.9940574363159407

Intrinsic single-mutation control:
- Spearman: 0.2525025092911592
- Pearson: 0.2991273011097151

After low-order nuisance residualization:
- partial context Spearman: -0.028251566611507472
- partial single Spearman: 0.048849795187586476
- |context| - |single|: -0.020598228576079004

Permutation null:
- permutations: 1,000
- empirical two-sided p: 0.7202797202797203
- null mean: -0.0005902378401106802
- null SD: 0.07689197383413558
- null |rho| p95: 0.15285635325394029
- null |rho| p99: 0.19092228318368012

Primary quad gate:
- |partial rho| >= 0.15: FAIL
- permutation p <= 0.01: FAIL
- beats single-control |rho| by >= 0.05: FAIL

## Results — GB1 order 3

n = 192

Raw candidate-conditioned context score:
- Spearman: -0.23622263516262923
- Pearson: -0.1602520416567415

Evaluation-only affine fit:
- intercept: -0.8203019416850826
- slope: -0.5928478575431827
- RMSE: 1.5410937397232822
- target SD: 1.5612714110613937
- normalized RMSE: 0.9870761283431213

Intrinsic single-mutation control:
- Spearman: 0.13062216857011097
- Pearson: 0.06538540108329974

After low-order nuisance residualization:
- partial context Spearman: -0.006227450288907577
- partial single Spearman: 0.01652571684344736
- |context| - |single|: -0.010298266554539783

Permutation null:
- permutations: 1,000
- empirical two-sided p: 0.9310689310689311
- null mean: 0.0028277371619238813
- null SD: 0.07268493680971162
- null |rho| p95: 0.145645955972113
- null |rho| p99: 0.18410720098743996

## Results — GB1 order 4

n = 192

Raw candidate-conditioned context score:
- Spearman: -0.22074980332582808
- Pearson: -0.16530090124991553

Evaluation-only affine fit:
- intercept: -2.33366043496361
- slope: -1.3771340524258957
- RMSE: 3.0677461800729393
- target SD: 3.1105372784415346
- normalized RMSE: 0.9862431809883228

Intrinsic single-mutation control:
- Spearman: -0.019753886010362695
- Pearson: 0.009112457445335448

After low-order nuisance residualization:
- partial context Spearman: -0.03728501478447224
- partial single Spearman: 0.009202994872907795
- |context| - |single|: 0.02808201991156445

Permutation null:
- permutations: 1,000
- empirical two-sided p: 0.6113886113886113
- null mean: -0.0037916220193690148
- null SD: 0.07403765982043815
- null |rho| p95: 0.14539587458969697
- null |rho| p99: 0.1871617204242737

Primary quad gate:
- |partial rho| >= 0.15: FAIL
- permutation p <= 0.01: FAIL
- beats single-control |rho| by >= 0.05: FAIL

## Quad decision

Partial candidate-conditioned context Spearman:
- TrpB: -0.028251566611507472
- GB1: -0.03728501478447224

Same sign:
true

Both satisfy the preregistered weak-rejection condition:
true

Therefore:
PRETRAINED_SEQUENCE_CONTEXT_SCALAR_REJECTED

## Low-order residualization numerical diagnostics

The 16-column design (intercept + 15 frozen low-order features) had rank 13 in all four cells because several exact low-order aggregate features are algebraically dependent. This was handled by numpy least-squares and was not changed post hoc.

TrpB order 3:
- design rank: 13
- min singular value: 1.2475719766984784e-15
- max singular value: 165.45590847235724

TrpB order 4:
- design rank: 13
- min singular value: 2.0343770774255594e-15
- max singular value: 146.0468476826814

GB1 order 3:
- design rank: 13
- min singular value: 9.553969147218e-17
- max singular value: 86.29820551824818

GB1 order 4:
- design rank: 13
- min singular value: 4.327303957477843e-16
- max singular value: 111.61373077548068

## Immutable artifacts

Prediction-visible artifact:
- ID: 11169055105
- digest: sha256:2b823698b53dec82e4a05f0c5de99b8e7203e005253f339f9564ce25c60b8b34

Frozen-score artifact:
- ID: 11168806614
- digest: sha256:3abcfe0b0d279f038096d330c0dd9ad998ec64af61060326efaa0f3b1b89b25e

Sealed-reveal artifact:
- ID: 11168731241
- digest: sha256:b52d9091fc05d64439cb22f36fd7b08f9acae9c6f261ef2c178c3b46498da3df

Complete evidence bundle:
- ID: 11168668376
- digest: sha256:66bd9258f4879491d61cc45713d6b54a1f609a874fcc85e19bd64b793b67d9b4

## Canonical file hashes

Prediction-visible:
- SELECTED_CANDIDATES.csv: 2749db7b5f25ded4007b8c26cce59cf09083eff888af797294ec7f159f8ebfaf
- LOW_ORDER_DIAGNOSTICS.json: fe425036e9d6381430e5b089b8340f80e74a200f467170b3a345a494039f5f86
- SOURCE_MANIFEST.json: 8a99d8b0bd620e38943eefab6fdfd1473b24e522c2e187ba9e02bade018b5d08

Frozen scores:
- ESM_CANDIDATE_SCORES.csv: 18a41f56c6f835e4f4f6fce44dac85fe8a4532120ebbc9d3815b96e2530e3a3b
- ESM_MUTATION_SCORES.csv: 46d1a082044188ceea09d9b45458d13007c885586e312c1b7e9c763810c6a624
- MODEL_MANIFEST.json: 08ecfecdb3b2ced90394980f3877306dd4e2e0c2a3a9020e4597e38873d12d99
- RUNTIME_MANIFEST.json: 1771995bef5bdbc7bcb7fdd13e76b324f82aa111e4ee418fd1b2969852173390
- SCORING_DIAGNOSTICS.json: 11506d3459c6c87d8f1dc382f0e096013dc7df27104e8be54f28f511bc542826

Sealed reveal:
- REVEAL_TARGETS.csv: 43a97dcda66bc9ec8a708dd4babcd3b290103137e82a5c527006193e516ae26d

Evaluation:
- CELL_METRICS.json: 42a6dbe56e9aab97d511057ce6bab41a465acddd4ae235dbcd3696e07c2f5960
- DECISION.json: c09b18a794b8af1140004fdf171d2cec73fc019ef1f34a6825d230848adcdced
- EVALUATION_DIAGNOSTICS.json: 12521cc1bd2001fb95a10faf3e6dcb666a88bad8c6a22044a3964ead15f2150b
- JOINED_EVIDENCE.csv: 499fc6422b89114223365273dca83c7b0772f123ff4b4adae4cf8053b06a122d
- PERMUTATION_NULL.csv: 859f24a8eb003fe2c6e5c8d94bdb4b3c2f87e1dc5574902c1d9da50e7bcfe4a2

## Interpretation boundary

The raw ESM context score showed moderate negative associations in some cells, especially GB1, but those associations collapsed after the frozen low-order nuisance features were removed.

On the primary quad question:
- TrpB partial rho = -0.0283, p = 0.7203
- GB1 partial rho = -0.0373, p = 0.6114

These values are compatible with the permutation null and do not beat the intrinsic single-mutation control by the preregistered margin.

Therefore this scalar candidate-conditioned ESM preference-shift score does not provide evidence of the missing transferable higher-order information.

Do not:
- build an architecture from this scalar score;
- tune ESM model size using these same results;
- tune score weights or thresholds;
- reinterpret the raw GB1 correlations as a PASS;
- consume NucB.

The broader question of whether a richer frozen sequence representation contains useful information remains unresolved and requires a separately preregistered hypothesis, not a repair of this score.
