# 2A.1 — Baselines

This stage compares acquisition controls using the leak-safe 2A.0 campaign simulator.

Read `PREREGISTRATION.md` before interpreting results.

## Baselines

- deterministic random;
- greedy current V8.3;
- exploration-only;
- historical fixed 50/50;
- static initial V8.3.

## Default run

The default RhlA development tournament uses:

- 5% common bootstrap;
- 20% final measurement budget;
- batch size 16.

The full landscape truth is available only to the virtual assay and evaluator. Policies never receive unrevealed labels.

## Outputs

Each baseline writes:

- campaign audit manifests;
- `DISCOVERY_CURVE.csv`;
- `BASELINE_SUMMARY.json`.

The tournament writes:

- `BASELINE_COMPARISON.csv`;
- `TOURNAMENT_MANIFEST.json`.

Results belong only under `results/`.
