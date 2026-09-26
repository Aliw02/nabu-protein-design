# NABU V9 TrpB Two-to-Many Blind External Validation — Result

Status: COMPLETE
Decision: TRPB_BLIND_EXTERNAL_TRANSFER_FAIL
Phase 3: CLOSED
Post-reveal tuning on TrpB test targets: FORBIDDEN

## Canonical execution

GitHub Actions run:
- 36243991521

Final result artifact:
- ID: 10907556059
- digest: sha256:ea4fe1342c2d561d8574e21da1854aafa0549f2abebf132c1ada557867972a23

Frozen pre-reveal prediction artifact:
- ID: 10907551005
- digest: sha256:3a3741bb35d4ba506d71f544607776ab8bb61ae02ffedbcebc6fc244125b2dcb

## Blind integrity

- Source MD5 matched: a611408d2db06907a023ddc8a1d94c12
- Source SHA256: 2676d06452f4f5dbd52cb5c8cad17d34c94ac8307bda581b527be208a3925ba5
- Fit: 8,633
- Validation: 2,158
- Test: 217,507
- Validation targets used for fit: false
- Test targets visible during prediction: false
- Freeze completed before reveal: true
- Frozen prediction semantic SHA256:
  ebd11befb2a971b2b984f3fd48db08afe746dae7d6e5af7e357a28f4fd17f1d0
- Consensus rank SHA256:
  e93486e2c284758ecd439ac138994ba8b5ff303292ace7a188b22ef7de4a1db3

## Train-only diagnostic

- clean-double OOF Spearman: 0.8784987991
- clean-double OOF RMSE: 0.1657255258
- all OOF Spearman: 0.8888781223
- selected Ridge alpha: 1.0

This strong low-order diagnostic did not translate into superior global
higher-order ranking on the blind test.

## Full blind test

| Arm | Spearman | NDCG | Top-1% hits | Regret |
|---|---:|---:|---:|---:|
| B2_MAIN_ONLY | **0.4706975858** | 0.9932455220 | 926 / 2176 | 0.0 |
| V9_CONTEXTUAL_ESM_RESIDUAL | 0.2439605939 | 0.9924142745 | **944 / 2176** | 0.0 |
| V9_CONSENSUS_B2_CONTEXT | 0.4669718433 | **0.9932510959** | 932 / 2176 | 0.0 |

Published reference points preregistered before reveal:

| Reference | Spearman | NDCG |
|---|---:|---:|
| FLIP2 Ridge one-hot | 0.437 | 0.994 |
| FLIP2 CARP-640M supervised | 0.509 | 0.996 |

Observations:
- B2 alone exceeds the published Ridge Spearman reference.
- Consensus also exceeds the published Ridge Spearman reference.
- Consensus does not beat B2 on full-test Spearman.
- Consensus NDCG is slightly above B2.
- Context alone improves Top-1% hits over B2, despite much weaker global rank correlation.

## Complete-clean-main region

Rows: 105,813

| Arm | Spearman | NDCG | Top-1% hits | Regret |
|---|---:|---:|---:|---:|
| B2_MAIN_ONLY | **0.4831878855** | 0.9893707022 | 315 / 1059 | 0.0 |
| V9_CONTEXTUAL_ESM_RESIDUAL | 0.0483494929 | 0.9881219499 | **389 / 1059** | 0.0 |
| V9_CONSENSUS_B2_CONTEXT | 0.4755661460 | **0.9895965470** | 331 / 1059 | 0.0 |

The contextual residual strongly changes elite selection but damages global
ordering within the region where it is active.

## Mutation-order split

### Order 3
- B2 Spearman: 0.6751217397
- Context Spearman: 0.3845811908
- Consensus Spearman: 0.6606472967

### Order 4
- B2 Spearman: 0.1778461371
- Context Spearman: 0.1364628942
- Consensus Spearman: 0.1777623376

Consensus remains positive at both orders, but does not exceed B2.

## Preregistered gate

PASS:
- Top-1% hits >= B2
- regret <= B2
- positive order-3 Spearman
- positive order-4 Spearman
- frozen hashes present
- source/split contracts pass
- validation targets unused

FAIL:
- Consensus full-test Spearman > B2
- Consensus clean-region Spearman > B2

Therefore:

    TRPB_BLIND_EXTERNAL_TRANSFER_FAIL

The published Ridge threshold alone is exceeded by both B2 and Consensus, but
the published-benchmark gate was defined to require the internal transfer gate
first, so it also fails formally.

## Scientific interpretation

The blind result does not support retaining the IRED-developed Pareto Consensus
as a generally superior higher-order transfer architecture.

However, the result reveals three useful facts:

1. The simpler B2 baseline generalizes strongly to TrpB and exceeds the
   preregistered published Ridge Spearman reference.
2. The frozen sequence-context residual contains substantial elite-selection
   signal (944 Top-1% hits vs 926 for B2) even though its global ranking is much
   worse.
3. The fixed Pareto consensus rule preserves most of B2's global ranking and
   slightly improves NDCG/Top-1 hits, but does not improve Spearman over B2.

TrpB test targets are now development-revealed and must not be used to tune a
replacement architecture for a future claim. Any new architecture requires a
fresh untouched external benchmark for authoritative validation.
