# NABU V10.0 Interaction-Shape Prior — Development Result

Date: 2026-09-28
Branch: exp/nabu-v10-interaction-shape-dev
Run: 36418960080
Artifact: 10968927659
Artifact digest: sha256:ff09fd290a2da7352787ab4764d5af4dee234535506981999f8f6edc62ed36fe

Decision:
V10_DEV_REJECT

Architecture:
B2_PLUS_EXTERNAL_INTERACTION_SHAPE_RANK_CORRECTION

Kill rule:
APPLIES

NucB consumed: false
Phase 3 opened: false

## GB1 -> TrpB

Support:
- V10-supported: 39,511 / 217,507 = 18.1654%
- unsupported rows remain exactly B2

Whole target:
- B2 Spearman: 0.4706975858
- V10 Spearman: 0.4435388023
- O2-only control: 0.4692246972
- shuffled control: 0.4704635451

- B2 NDCG: 0.9932455220
- V10 NDCG: 0.9909408910

- B2 Top-1% hits: 926
- V10 Top-1% hits: 348

Order 3:
- B2 Spearman: 0.6751217397
- V10 Spearman: 0.6031789761

Order 4:
- B2 Spearman: 0.1778461371
- V10 Spearman: 0.1739216711

Supported region:
- B2 Spearman: 0.5630724957
- V10 Spearman: 0.5546910059

Gate failures:
- whole Spearman did not beat B2
- NDCG fell below B2
- Top-1% hits fell below B2
- order-3 Spearman fell below B2
- order-4 Spearman fell below B2
- V10 did not beat O2-only control
- V10 did not beat shuffled control

Deterministic replay: exact
Unsupported fallback: exact B2

## TrpB -> GB1

Support:
- V10-supported: 137,321 / 147,193 = 93.2932%
- unsupported rows remain exactly B2

Whole target:
- B2 Spearman: 0.3311304363
- V10 Spearman: 0.3322407573
- O2-only control: 0.3374907865
- shuffled control: 0.3298734597

- B2 NDCG: 0.7742795543
- V10 NDCG: 0.8155348683

- B2 Top-1% hits: 240
- V10 Top-1% hits: 336

Order 3:
- B2 Spearman: 0.5151622384
- V10 Spearman: 0.4794760202

Order 4:
- B2 Spearman: 0.2958106709
- V10 Spearman: 0.2949219489

Supported region:
- B2 Spearman: 0.3407455948
- V10 Spearman: 0.3504776168

Gate failures:
- order-3 Spearman fell below B2
- order-4 Spearman fell below B2
- V10 did not beat O2-only control

Deterministic replay: exact
Unsupported fallback: exact B2

## Frozen interpretation

The final ablation established that transferable interaction-shape information exists beyond scalar O2 cancellation.

V10.0 tested one explicit use of that information:
external interaction-shape prediction of B2 global rank error.

That architecture does NOT transfer robustly.

It damages TrpB substantially and, although it improves GB1 aggregate NDCG/Top-1% and slightly improves global Spearman, it degrades both mutation-order Spearman metrics and loses to the O2-only correction control.

Therefore:
- V10.0 is rejected;
- do not patch its blend, weight, thresholds, support rule, or order handling using these revealed results;
- do not create V10.1 as a repair of this same architecture;
- do not consume NucB;
- Phase 3 remains closed.

The interaction-shape signal remains a measured property of the revealed landscapes, but this rank-correction architecture is not a valid general solution.
