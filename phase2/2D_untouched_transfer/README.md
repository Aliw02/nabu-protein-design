# 2D — Untouched External Benchmark

This stage is quarantined until 2C is frozen.

Read first:

- `../FINAL_VALIDATION_PROTOCOL.md`
- `../BENCHMARK_LEDGER.json`
- `../RESULTS_REGISTRY.md`

Reserved blind evaluations:

1. **FLIP AAV/Random**
   - direct active-optimization head-to-head;
   - NABU vs greedy vs UCB vs Thompson Sampling vs random;
   - primary metric: AUC of Top-100 recovery versus measurements.

2. **FLIP2 IRED two-to-many**
   - low-order to higher-order mutation transfer;
   - comparison against published FLIP2 prediction baselines where metrics are directly comparable.

Do not place blind dataset labels or downloaded blind data in this directory before the 2C freeze.

Final outputs belong under:

`results/final_blind/`

All failed runs and audited reruns must be preserved.
