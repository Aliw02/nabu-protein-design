# NABU V8.1 GB1 Development Replication Protocol

Status: DEVELOPMENT REPLICATION — NOT A CANONICAL FREEZE.

## Purpose

V8.1 showed that the adaptive higher-order signal on PHOT survives strict cross-fitting. This replication asks whether the exact same V8.1 architecture and hyperparameters transfer to a second combinatorial protein landscape without changing the model.

The model runner is not modified for this replication.

## Dataset

Source:
- Repository: rcronk/mutate
- Commit: 0861a3f1f6c58a9597fde689eb7f9a07b38e2638
- Path: ridge/data/gb1_wu2016_landscape.csv
- Columns: variant, fitness
- Landscape: four-site GB1 Wu 2016 combinatorial fitness landscape

This dataset has been seen in prior project experiments, so this run is replication/development evidence only, not fresh blind confirmation.

## Representation adapter

The source encodes each genotype as four amino-acid letters. The V8.1 engine expects mutation-like categorical tokens.

The adapter maps a genotype ABCD to:

X1A:X2B:X3C:X4D

where X is only a neutral fixed prefix. This does not inject amino-acid properties, structural information, fitness values, or biological priors. It simply represents each position-state as a categorical token.

All variants therefore have four categorical states and pass the existing 3-5 token scope unchanged.

## Frozen architecture

Use the exact runner already committed at:

experiments/v8_1_crossfit_higher_order_dev/run_v8_1_crossfit_dev.py

No changes to:
- 5-fold cross-fitting
- outer visible/hidden FNV split
- additive main shrinkage
- raw-pair layer
- triplet support >= 2 and n/(n+2)
- quartet support >= 2 and n/(n+4)
- adaptive support coverage
- global permutation control
- mutation-count-stratified permutation control
- bootstrap configuration

## What we want to know

Replication is encouraging if:
- B5 Cross-Fit Adaptive Higher-Order remains above B3 Raw Pair on hidden Spearman.
- The direction is positive across the four-state landscape rather than only a tiny candidate subset.
- Negative controls collapse.
- Candidate ranking changes meaningfully.
- Bootstrap delta is positive, while recognizing that row bootstrap is not a substitute for an independent assay.

If the higher-order gain disappears or reverses, do not retune this branch. Record it as evidence that PHOT-specific higher-order structure did not transfer under the frozen V8.1 architecture.

## Breakpoint rule

Do not create a canonical architecture freeze until PHOT V8.1 and this second-landscape replication have both been reviewed together.

If the mechanism transfers, the V8.1 runner becomes a breakpoint candidate and the next stage is a sealed fresh benchmark on an unrevealed landscape.
