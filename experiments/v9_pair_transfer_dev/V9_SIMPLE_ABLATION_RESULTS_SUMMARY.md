# NABU V9 Simple Architecture Ablation — Result

Status: COMPLETE
Ablation winner: V9_CONSENSUS_B2_CONTEXT
Phase 3: CLOSED
Scientific claim: false

## Canonical execution

GitHub Actions run:
- 36238279578

Workflow head:
- eeabc2ef4811abf64776fff06eda84925aafb3ba

Artifact:
- 10904602661
- sha256:5baccb85e989ba5e368150ac846a024e4fb701e1b69841355bb5dc9de03090b0

Pinned IRED source:
- aa45a2f85fb1af87b6b0e86397b3f8f292a061dc574be2536673b5bd490b0e74

All preregistered contracts passed.

## Ablation

Two simple alternatives were tested:

A. V9_CONSENSUS_B2_CONTEXT
- no fitted weights;
- no scalar blending;
- deterministic Pareto layers over B2 and Context rankings;
- arithmetic mean fractional rank only as within-layer tie-break.

B. V9_CONTACT_AWARE_CONTEXT
- frozen ESM-2 sequence context;
- one additional 320-dimensional contact-neighborhood channel;
- same fit/train-only RidgeCV residual readout;
- no distance cutoff or hand-authored biochemical rules.

## Full IRED development test

| Arm | Spearman | NDCG | Top-1% hits | Normalized regret |
|---|---:|---:|---:|---:|
| B2_MAIN_ONLY | 0.1699713866 | 0.9560896339 | 3 / 42 | 0.0647016495 |
| V9_CONTEXTUAL_ESM_RESIDUAL | 0.1734515583 | 0.9564649249 | 1 / 42 | 0.0647016495 |
| V9_CONSENSUS_B2_CONTEXT | **0.1769446973** | **0.9569780968** | **3 / 42** | **0.0647016495** |
| V9_CONTACT_AWARE_CONTEXT | 0.1725728223 | 0.9565198335 | 1 / 42 | 0.0647016495 |

Consensus is the only arm that simultaneously:
- improves full-test Spearman over B2;
- preserves B2's 3/42 Top-1% hits;
- preserves B2's normalized regret;
- improves NDCG over both B2 and Context.

## Complete-clean-main region

Rows: 1,871

| Arm | Spearman | NDCG | Top-1% hits | Normalized regret |
|---|---:|---:|---:|---:|
| B2_MAIN_ONLY | 0.2114864341 | 0.9480111432 | 3 / 19 | 0.0593089448 |
| V9_CONTEXTUAL_ESM_RESIDUAL | 0.2572253615 | 0.9504729130 | 1 / 19 | 0.0593089448 |
| V9_CONSENSUS_B2_CONTEXT | **0.2271452099** | 0.9495277368 | **3 / 19** | **0.0593089448** |
| V9_CONTACT_AWARE_CONTEXT | 0.2554801603 | **0.9506008696** | 1 / 19 | 0.0593089448 |

Context alone remains strongest for clean-region Spearman, but loses elite overlap.
Consensus deliberately trades part of that rank gain to preserve B2's elite behavior while still beating B2 on clean-region Spearman.

## Mutation-order stability

Consensus Spearman:
- order 3: 0.2380901717
- order 4: 0.1787552958
- order 5: 0.0306011353
- order 6: 0.0475909197
- order 7: 0.0753319714

The preregistered positive Spearman requirement for orders 3, 4, and 5 passes.

## Contact-aware residual diagnostic

Frozen contact-aware readout:
- selected Ridge alpha: 1000.0

Fit/train OOF:

Base Context:
- Spearman: 0.1988308964
- RMSE: 1.0077471657

Contact-aware:
- Spearman: 0.1966765167
- RMSE: 1.0137536525

The added contact-neighborhood channel does not improve the residual-learning diagnostic and does not repair Top-1% selection.

## Preregistered gates

### V9_CONSENSUS_B2_CONTEXT

PASS:
- full Spearman >= B2
- clean-region Spearman >= B2
- Top-1% hits >= B2
- regret <= B2
- positive Spearman at orders 3/4/5
- deterministic replay

Arm result: PASS

### V9_CONTACT_AWARE_CONTEXT

PASS:
- full Spearman >= B2
- clean-region Spearman >= B2
- regret <= B2
- positive Spearman at orders 3/4/5
- deterministic replay
- residual-training OOF Spearman > 0

FAIL:
- Top-1% hits >= B2

Arm result: FAIL

## Predeclared ablation decision

    V9_CONSENSUS_B2_CONTEXT

is the ablation winner.

## Architectural interpretation

The simplest successful solution is not to force B2 and Context into one
fitness scalar.

B2 retains useful elite-selection information while the Context model adds
higher-order ranking information. Deterministic Pareto consensus preserves both
signals without a learned weight or post-hoc scalar.

The tested contact-aware channel is unnecessary in its present form.

This is still post-Phase-2 development evidence on already-revealed IRED.
It does not repair the frozen Phase-2 verdict and does not open Phase 3.
The winning V9 consensus architecture requires a new untouched external
benchmark for scientific validation.
