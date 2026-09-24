# NABU V8.2 Development Tournament

This branch is a development descendant of the frozen V8.1 breakpoint:

- Base branch: `freeze/nabu-v8-1-crossfit-higher-order-breakpoint`
- Base commit: `bce39ee5cbcb05259375a2c1f212a49f69834bdb`
- Frozen branch is not modified.

## Purpose

V8.1 showed an objective split:

- PHOT: higher-order improved whole-landscape ranking and elite selection.
- GB1: triplet higher-order slightly harmed whole-landscape Spearman but improved elite Top-K selection.

V8.2 is intentionally a **revealed development tournament**, not a prospective scientific confirmation. Its purpose is to diagnose how much higher-order signal should be allowed to change the pairwise ranking and where.

## Shared invariant

Every arm uses the same V8.1 cross-fitted hierarchy. The tournament changes only the final decision/scoring policy, so failures can be attributed to the policy rather than a different training split or a different memory fit.

## Tournament directions

### Baselines

- `B3_RAW_PAIR`
- `B4_CROSSFIT_TRIPLET`
- `B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER`

### A. Higher-order dosage

- `T_QUARTER_HIGHER_ORDER`
- `T_HALF_HIGHER_ORDER`
- `T_THREE_QUARTER_HIGHER_ORDER`

Question: is the GB1 global-rank loss simply excessive higher-order correction?

### B. Stronger confidence shrinkage

- `T_CONFIDENCE_SQUARED`

Question: should lower-support higher-order landscapes be shrunk more aggressively?

### C. Elite-only objective split

- `T_ELITE_10PCT_HIGHER_ORDER`
- `T_ELITE_20PCT_HIGHER_ORDER`
- `T_ELITE_30PCT_HIGHER_ORDER`

Question: can higher-order be used where protein design actually cares most while leaving most global pairwise ordering untouched?

### D. Scale-free pair / higher-order rank fusion

- `T_RANK_FUSION_50_50`
- `T_RANK_FUSION_PAIR_75`

Question: is score-scale mismatch responsible for some of the cross-landscape instability?

### E. Visible-OOF objective gate

- `T_OOF_OBJECTIVE_GATE`

This is the main adaptive hypothesis.

The gate uses only visible OOF behavior:

- if B4 OOF Spearman > B3 OOF Spearman, use higher-order globally;
- otherwise use higher-order only inside the B3 top 20% region.

The hidden truth is not used to choose the gate mode.

## Diagnostics

For every arm and landscape the runner stores:

- whole-landscape Spearman;
- Top-1 true percentile;
- Top-5 / Top-10 / Top-50 mean true percentile;
- Top-10 and Top-50 true Top-1% hit counts;
- Top-50 enrichment;
- Top-50 mean true fitness;
- normalized Top-1 regret;
- Top-5 / Top-10 / Top-50 candidate replacement analysis versus B3;
- mean and p95 rank shift versus B3;
- mean and p95 score perturbation versus B3;
- diagnostic tags for global-rank gain/harm and elite gain/no-gain.

The replacement analysis is especially important: it reports whether candidates promoted into Top-K are actually stronger than the B3 candidates they displaced.

## Run

From repository root:

```powershell
py -m pip install pandas numpy scipy

py experiments/v8_1_gb1_replication_dev/prepare_gb1.py `
  gb1_wu2016_landscape.csv `
  --out GB1_Wu2016_NABU_V8_1.csv

py experiments/v8_2_objective_tournament/run_v8_2_tournament.py `
  PHOT_CHLRE_Chen_2023.csv `
  GB1_Wu2016_NABU_V8_1.csv `
  --out nabu_v8_2_tournament_results
```

## Outputs

- `nabu_v8_2_tournament_results/V8_2_TOURNAMENT_SUMMARY.json`
- `nabu_v8_2_tournament_results/V8_2_TOURNAMENT_ALL_ARMS.csv`
- one subdirectory per dataset containing:
  - `V8_2_DATASET_RESULTS.json`
  - `V8_2_ARM_SUMMARY.csv`
  - `V8_2_CANDIDATE_SCORES.csv`

## Development rule

Do not freeze V8.2 merely because one arm wins one dataset. Use the tournament to identify a stable mechanism across PHOT and GB1, then add more already-revealed compatible landscapes if available. Only after the policy is stable should a future untouched benchmark be reserved for prospective confirmation.
