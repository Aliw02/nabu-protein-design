# Phase 2A.4 — Multi-Landscape Evaluation Preregistration

## Goal

Test whether the acquisition behavior developed on RhlA transfers across multiple non-RhlA protein fitness landscapes without per-landscape retuning.

## Landscapes

Pinned historical sources already used by NABU Phase 1:

- GB1 (Wu 2016);
- TrpB (Johnston 2024);
- PhoQ.

These are development/generalization landscapes, not the final untouched Phase-2D transfer set.

## Controllers frozen before evaluation

1. `deterministic_random` — sanity floor.
2. `historical_50_50` — fixed reference.
3. `relative_evidence_v2` — adaptive candidate frozen after 2A.3.

No controller rule may inspect the landscape name.

No controller parameter may be changed between landscapes.

## Common deterministic campaign pool

The source landscapes are much larger than the RhlA development pool.

To keep closed-loop computation comparable and reproducible, each landscape is reduced to a deterministic **identity-only** campaign subpool.

Rules:

- canonicalize candidate identities first;
- choose exactly 2048 candidates by a fixed FNV-1a identity hash;
- do not use fitness to construct the campaign pool;
- evaluate discovery relative to the truth distribution inside that frozen campaign pool.

Because the four-site landscapes use the same canonical genotype representation, identity-only hashing also maximizes protocol comparability.

## Shared measurement protocol

For every landscape and controller:

- campaign-pool size: 2048;
- bootstrap: 64 measured candidates;
- final budget: 256 measured candidates;
- frozen Phase-1 V8.3 core;
- same identity-only bootstrap logic;
- same oracle isolation;
- same metrics;
- deterministic tie-breaking.

Primary cadence:

- micro_8.

Sensitivity cadence:

- micro_16.

## Experiment matrix

Primary micro_8:

- deterministic_random;
- historical_50_50;
- relative_evidence_v2.

Sensitivity micro_16:

- historical_50_50;
- relative_evidence_v2.

Total: five conditions per landscape.

## Primary metrics

- normalized post-bootstrap discovery AUC;
- final best true fitness;
- best-so-far regret;
- cumulative Top-1% hits;
- first Top-1% measurement.

Diagnostics:

- scoreable coverage;
- visible B3 OOF Spearman;
- router mode;
- V2 exploration-pressure trajectory.

## Generalization interpretation

The main comparison is paired within each landscape:

`relative_evidence_v2 - historical_50_50`

at identical cadence and budget.

Random is a sanity floor, not a controller-development target.

Do not choose a new controller formula after seeing one landscape and then count the other landscapes as untouched evidence in the same run.

## Gate to 2B

2A.4 is sufficient to move toward 2B when we can state, across multiple landscapes:

- whether closed-loop acquisition beats or meaningfully differs from random;
- whether historical 50/50 is stable;
- whether Relative-Evidence V2 generalizes, is neutral, or degrades;
- which cadence/controller pair is defensible as the Phase-2 acquisition reference.

A universal superiority claim is not required to proceed.
