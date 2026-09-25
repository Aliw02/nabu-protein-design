# Phase 2A.2 Micro-Batch Refit Ablation — Preregistration

## Question

Does refitting more frequently improve discovery when the acquisition policy itself is held fixed?

This stage tests **batch staleness only**.

## Fixed acquisition policy

Use the frozen Phase-2A.1 `historical_50_50` policy without changing:

- structural exploration definition;
- exploitation ranking;
- 0.5 / 0.5 rank fusion;
- tie-breaking;
- frozen V8.3 core.

No new adaptive controller is allowed in this stage.

## Development landscape

RhlA training pool as a virtual assay oracle.

No RhlA hidden-triple test truth is used for acquisition.

## Shared protocol

All cadence conditions use:

- identical identity-only 5% bootstrap;
- identical bootstrap IDs;
- identical candidate pool;
- identical oracle;
- identical final 20% measurement budget;
- identical acquisition policy;
- identical metrics.

Pool size: 1608.

Default counts:

- 5% seed: 81 measurements;
- 10% checkpoint: 161 measurements;
- 20% final budget: 322 measurements.

## Cadence conditions

1. **historical_jumps**
   - refit after reaching 10%;
   - then refit after reaching 20%;
   - additions: 80, then 161.

2. **micro_32**
   - refit every 32 revealed measurements.

3. **micro_16**
   - refit every 16 revealed measurements.

4. **micro_8**
   - refit every 8 revealed measurements.

5. **micro_4**
   - refit every 4 revealed measurements.

The final partial batch, if any, is truncated so every condition stops at exactly 322 total measurements.

## Primary metrics

- normalized post-bootstrap discovery AUC;
- final best true fitness;
- final regret;
- cumulative Top-1% discoveries;
- first Top-1% measurement checkpoint.

Diagnostics:

- final scoreable coverage;
- final visible B3 OOF Spearman;
- router trajectory;
- number of refits.

## Interpretation boundary

This experiment can show whether refit cadence matters on RhlA under the fixed historical 50/50 policy.

It cannot establish a universal optimal batch size from one landscape.

Do not tune the Phase-2A.3 adaptive controller from one preferred RhlA batch size and report that as untouched evidence.
