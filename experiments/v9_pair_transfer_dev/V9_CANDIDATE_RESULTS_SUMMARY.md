# NABU V9 Clean Calibrated Pair Expansion — Result

Status: COMPLETE
Decision: V9_CLEAN_CALIBRATED_PAIR_EXPANSION_REJECT
Phase 3: CLOSED
Scientific claim: false

## Execution

GitHub Actions run:
- 36235511196

Workflow head:
- 5e270e470748d34d28a098a45d27d95a50e29847

Artifact:
- 10904400780
- sha256:b9c0e508e500a956518c1cccd94275f86bfc0e5ac359218b7b64cfc4af48fed6

All contract tests and the workflow completed successfully.

## Core comparison on IRED development test

| Arm | Spearman | Top-1% hits | Normalized regret |
|---|---:|---:|---:|
| Phase-2 V8.3 abstention | -0.0217706 | 0 / 42 | 0.2315224 |
| B2 main-only | 0.1699714 | 3 / 42 | 0.0647016 |
| V9.0 raw pair transfer | 0.0951920 | 0 / 42 | 0.2157016 |
| V9.1 degree-normalized transfer | 0.1370972 | 3 / 42 | 0.0647016 |
| V9 clean calibrated pair expansion | 0.1551567 | 0 / 42 | 0.2103395 |

The clean calibrated V9 candidate improves over V9.1 in full-test Spearman but
does not beat B2 and loses the elite-selection behavior.

## Coverage

- test rows: 4,178
- B2 all-main supported: 3,648
- complete clean-single support: 1,871
- B2 fallback rows: 2,307
- rows with exact clean pairs: 62
- rows with calibrated unseen-pair predictions: 1,861
- exact clean pair contributions: 63
- calibrated unseen-pair contributions: 8,637
- unresolved pair contributions: 17,559

## Clean-supported region

On the same 1,871 rows with complete clean-single support:

- B2 Spearman: 0.2114864
- V9.1 Spearman: 0.1669960
- V9 clean calibrated expansion Spearman: 0.1789433

Thus the full-test rejection is not caused only by the B2 fallback rows.

## Mutation-count behavior

V9 clean calibrated expansion Spearman:
- order 3: 0.2198943
- order 4: 0.1641791
- order 5: 0.0155791
- order 6: 0.0396371
- order 7: 0.0409712
- order 8: 0.0667260
- order 9: 0.4303030 (n=10)
- order 10: -0.2581989 (n=4)

The preregistered requirement of positive Spearman at orders 3/4/5 passed, but
performance decays sharply by order 5 and elite selection is poor.

## Post-run mechanism diagnostic

This diagnostic does not change the preregistered decision.

For the 1,871 complete-clean-main rows:

    clean_main_baseline = V9_score - summed_pair_correction
    true_missing_residual = target - clean_main_baseline

The summed calibrated pair correction retains real information:

- pair-sum vs true missing residual Spearman: 0.4743684
- Pearson: 0.4834595

By mutation order on clean-supported rows:

| Order | n | pair-sum vs residual Spearman | pair-sum std | residual std |
|---|---:|---:|---:|---:|
| 3 | 1111 | 0.3943175 | 1.02647 | 1.64852 |
| 4 | 487 | 0.4547242 | 1.82214 | 1.81409 |
| 5 | 193 | 0.4264212 | 2.57570 | 2.00504 |
| 6 | 53 | 0.5400742 | 3.32948 | 2.73504 |
| 7 | 16 | 0.3382353 | 5.39029 | 2.59202 |

Therefore the clean pair predictor is not devoid of transferable signal.
The failure occurs when pair effects are composed as a simple second-order
additive sum in a higher-order mutational background. Pair-sum scale grows with
order and the combined ranking/elite selection worsens despite the pair-sum
retaining correlation with the missing residual.

## Preregistered decision checks

PASS:
- full Spearman beats V9.1
- positive Spearman at mutation orders 3, 4 and 5
- prediction scale remains within 3x target standard deviation
- deterministic rank replay

FAIL:
- full Spearman does not beat B2
- Top-1% hits do not match B2
- normalized regret is worse than B2
- clean-supported-region Spearman does not beat B2

Decision:

    V9_CLEAN_CALIBRATED_PAIR_EXPANSION_REJECT

## Architectural conclusion

The clean-main and calibrated clean-pair discoveries are retained as useful
diagnostic components, but identity-only pair effects cannot simply be summed
to produce a general higher-order score.

Per the preregistered boundary, stop scalar, normalization, and pair-sum tuning
for this representation family.

The next V9 architecture must represent interaction context directly using a
richer relational representation (sequence-context and/or structure-aware
features) rather than treating higher-order fitness as a context-free sum of
identity-only pair epistasis.

IRED remains development-only. Phase 3 remains closed.
