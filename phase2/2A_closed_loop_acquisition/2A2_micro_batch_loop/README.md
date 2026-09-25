# 2A.2 — Micro-Batch Refit Ablation

This stage isolates **batch staleness**.

The acquisition rule is held fixed as the preregistered historical 50/50 policy from 2A.1.

Only the reveal/refit cadence changes.

## Conditions

- historical 5% → 10% → 20% jumps;
- fixed batches of 32;
- fixed batches of 16;
- fixed batches of 8;
- fixed batches of 4.

Every condition starts from the same 5% bootstrap and ends at the same 20% total measurement budget.

## Outputs

Each cadence writes:

- campaign manifests;
- discovery curve;
- cadence summary.

The experiment writes:

- `MICRO_BATCH_COMPARISON.csv`;
- `MICRO_BATCH_MANIFEST.json`.

Read `PREREGISTRATION.md` before interpreting results.
