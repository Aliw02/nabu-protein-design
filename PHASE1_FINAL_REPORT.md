# NABU Phase-1 Final Report

## Final status

`PHASE1_FINAL_GATE_CLOSED_WITH_LIMITATION`

Phase 1 is closed. The limitation applies to the experimental Active Acquisition controller, not to the frozen V8.3 core.

Final closure reference:

`freeze/nabu-phase1-final-cr9114-closed`  
Commit: `690494691c20edc3c9b7a6bba9a5156913173c58`

## Phase-1 objective

Phase 1 asked whether a discrete, CPU-oriented combinatorial model could learn useful mutation interaction structure from measured variants and generalize to unseen combinations without requiring protein structure coordinates or a large pretrained neural model.

The required capabilities progressively became:

- recover additive signal,
- model pairwise epistasis,
- model supported triplet/quartet residual structure,
- avoid same-row higher-order leakage,
- decide when higher-order terms help,
- protect a strong pairwise backbone when they do not,
- rank hidden combinatorial variants,
- operate under limited measurement budgets,
- explore whether active experiment selection can reduce measurements.

## Final frozen core

NABU V8.3 Dual-Objective Router.

Core freeze:

`freeze/nabu-v8-3-dual-objective-router-validated`  
`c8afdcd7698231d95a39ef13a3fe22b6e3f507c5`

The core contains:

- B2 additive main memory,
- B3 pairwise residual memory,
- B4 cross-fitted triplet residual memory,
- B5 cross-fitted quartet residual memory,
- visible OOF global gate,
- visible OOF elite gate,
- B3 protection fallback.

## Evidence chain

### Multi-landscape development

PHOT and TrpB showed global higher-order benefit.

GB1 showed that higher-order could be useful for elite discovery while harming global rank, motivating elite-only reranking.

PhoQ showed the need for an explicit protection mode.

### CreiLOV

A sealed low-order-to-4/5-mutant extrapolation test.

Result: strong B3 generalization, no supported higher-order memory, correct router abstention.

### eqFP611

Fresh sealed higher-order test.

Result: supported triplet/quartet memory and strict hidden-set gains on 4/5 primary endpoints.

### RhlA

Passive sample-efficiency benchmark established that the frozen practical threshold first appeared at 40%.

Retrospective Active Acquisition development later reached the same threshold at 20%, suggesting 2× measurement efficiency.

### CR9114-H1 final gate

Fresh one-shot final validation with natural 4/5-mutation hidden candidates.

Qualification before reveal:

- 546 passive triplet entries
- 1,698 passive quartet entries
- 545 active-20 triplet entries
- 1,024 active-20 quartet entries
- 1,673 fixed hidden variants
- 100% hidden scoreability coverage

Core result:

- B3 Spearman: 0.9271
- V8.3 Spearman: 0.9387
- 4/5 strict primary gains
- 5/5 no regression
- Core gate: PASS

Active result:

- Active 10% reached Passive 40% practical threshold
- 75% measurement reduction
- 4× efficiency factor
- Active 20% failed the preregistered target
- Active 40% was also not monotonic

Because Active 20% was the preregistered target, the overall formal result remains closed with limitation.

## What Phase 1 supports

Phase 1 supports the claim that the frozen NABU V8.3 core can:

- learn useful pairwise interaction structure,
- exploit supported higher-order residual structure,
- detect when global higher-order use is justified,
- restrict higher-order use to elite reranking when appropriate,
- suppress higher-order use when unsupported,
- generalize to hidden higher-order combinatorial variants in published benchmark landscapes.

## What Phase 1 does not support

Phase 1 does not establish:

- de novo full-protein sequence design,
- direct structure prediction,
- structure-conditioned protein generation,
- monotonic universal performance of the current Active Acquisition controller,
- wet-lab confirmation of a new NABU-generated protein outside the benchmark datasets.

## Closure rule

No Phase-1 dataset should be used for further post-reveal tuning of the frozen V8.3 core.

Phase 2 must treat the Phase-1 core as a frozen baseline.
