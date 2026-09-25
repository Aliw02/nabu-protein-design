# NABU Phase 2 Final Validation Protocol

## Purpose

Prevent development-set overfitting and ensure that the final Phase-2 scientific verdict is based on datasets that were not used to design, tune, debug, or select the NABU Phase-2 architecture.

No dataset used in Phase 1, 2A, 2B, or 2C development may count toward the final blind PASS.

## Dataset classes

### Development / diagnostic only

The following datasets are permanently excluded from the final blind PASS calculation because they have already influenced development, debugging, model selection, or interpretation:

- RhlA
- GB1
- TrpB
- PhoQ
- PHOT
- eqFP611
- CreiLOV
- CR9114-H1

Results on these datasets remain scientifically useful, but they are non-blind evidence.

### Blind external benchmark A — AAV/Random

Reserved dataset:

- FLIP AAV sampled/random protein fitness landscape.
- This dataset has not been used in NABU development.

Published benchmark reference:

- Greenman, Amini & Yang, PLOS Computational Biology (2025), "Benchmarking uncertainty quantification for protein engineering".

Direct comparison category:

- retrospective Bayesian optimization / active protein fitness optimization.

Published protocol to reproduce:

- start from a random 10% sample of the available AAV/Random training pool;
- compare greedy, UCB, Thompson Sampling, and random acquisition;
- perform 5 acquisition stages with training-set sizes equally spaced on a log scale;
- repeat across 3 folds / initialization seeds;
- primary published metric: percentage of the global Top-100 fitness sequences recovered versus fraction of training data observed.

The released authors' implementation or an exact audited reproduction must be run beside NABU. Published figure values alone are not sufficient when a reproducible implementation is available.

### Blind external benchmark B — FLIP2 IRED two-to-many

Reserved dataset:

- FLIP2 Imine Reductase (IRED), two-to-many split.
- 17,143 measured variants.
- Train domain: variants with 0, 1, or 2 mutations.
- Test domain: higher-order mutants.
- This dataset is not used for Phase-2 development.

Purpose:

- external extrapolation test from low-order evidence to higher-order protein variants;
- directly relevant to NABU's current combinatorial-design regime.

Shared metrics:

- Spearman rank correlation;
- NDCG;
- Top-k enrichment;
- normalized regret;
- Top-1% recall where definable.

Published FLIP2 baselines are retained for context, but this is a fitness-inference transfer benchmark rather than the direct active-learning head-to-head.

## Blindness / quarantine rules

Before any AAV or IRED target labels are loaded:

1. Phase 2C architecture must be frozen.
2. Maturity-gate rule must be frozen.
3. Assembly rule and search budget must be frozen.
4. Acquisition controller(s) included in the final system must be frozen.
5. Batch cadence must be frozen.
6. Metrics and primary winner rule must be frozen.
7. Seeds / fold-generation procedure must be frozen.
8. Tie-breaking and stopping rules must be frozen.
9. Data adapters must be tested only on synthetic or development data.

The reserved blind datasets must not be used for:

- hyperparameter tuning;
- architecture selection;
- gate calibration;
- batch-size selection;
- acquisition-rule selection;
- choosing among controller versions;
- debugging scientific behavior.

After reveal, only implementation bugs that are demonstrably label-independent may be fixed. Every such change requires:

- a new version;
- an audit note;
- a preserved failed run;
- rerunning the complete benchmark suite.

No result may be silently overwritten.

## Direct head-to-head winner rule

Primary direct comparison: AAV/Random Bayesian-optimization benchmark.

All methods receive the same:

- initial fold;
- candidate pool;
- measurement budget;
- acquisition checkpoints;
- oracle labels;
- evaluation checkpoints.

Competitors:

- NABU Phase-2 frozen system;
- published greedy baseline;
- published UCB baseline;
- published Thompson Sampling baseline;
- random baseline.

Primary score:

- area under the Top-100-recovery-versus-measurement curve.

Secondary shared scores:

- final percentage of Top-100 recovered;
- best normalized fitness found;
- measurements to first Top-100 hit;
- cumulative Top-100 discoveries.

A head-to-head winner is reported from the preregistered primary score.

Do not redefine the winner using a secondary metric after observing outcomes.

Report:

- mean across folds;
- standard deviation;
- per-fold values;
- paired deltas versus every baseline.

If statistical uncertainty makes a difference indistinguishable, report a tie / unresolved comparison rather than forcing a winner.

## Phase-2 final scientific gate

Phase 2 may receive a final scientific PASS only after:

- 2C is frozen without using the blind datasets;
- AAV direct benchmark is executed exactly once under the frozen protocol, except audited label-independent bug reruns;
- IRED external higher-order transfer benchmark is executed under its frozen split;
- all runs and failures are preserved;
- the direct published-baseline comparison is reported without selective metric choice;
- limitations are explicitly retained.

The final Phase-2 verdict must distinguish:

- engineering integrity;
- active optimization performance;
- higher-order transfer;
- assembly quality;
- failure cases.

## Relationship to Phase 3

These final tests evaluate NABU's current combinatorial / closed-loop protein engineering category.

They are not a de-novo full-sequence generation benchmark.

Phase 3 de-novo design will require a separate benchmark suite.
