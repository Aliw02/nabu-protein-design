# NABU V8.1 Cross-Fitted Higher-Order Development Protocol

Status: DEVELOPMENT ONLY — NOT A CANONICAL BREAKPOINT.

## Why V8.1 exists

V8 Dev showed a promising gain from triplet/quartet residual memory over the raw-pair baseline on PHOT, but the higher-order memories were learned from residuals produced on the same visible rows used to fit lower-order layers.

V8.1 tests whether that gain survives when the residual targets are cross-fitted.

## Fixed development data

PHOT_CHLRE_Chen_2023 is development data. Its outcomes were already revealed in V7/V8, so this run is not a fresh blind confirmation.

## Deterministic split

The outer visible/hidden split remains unchanged:

- FNV1a32(canonical mutation id) mod 10
- buckets 0..6 = visible
- buckets 7..9 = hidden

Inside the visible set, V8.1 assigns five deterministic cross-fit folds with:

- FNV1a32(candidate_id + "|NABU_V8_1_CF|161") mod 5

No fitness labels are used to create folds.

## Cross-fitted hierarchy

### B2 / B3

For each visible fold:

- Fit additive main effects and raw-pair residual memory on the other four folds.
- Predict the held-out fold.
- Store the held-out B3 residual y - B3_oof.

The final triplet memory is learned from these OOF B3 residuals.

### B4

To produce an unbiased B4 residual target for every visible row, each outer fold is excluded completely.

For an outer fold f:

1. Use only the other four folds.
2. Within those four folds, generate inner OOF B3 residuals by excluding each inner fold g in turn.
3. Build triplet memory from those inner OOF B3 residuals.
4. Predict the untouched outer fold f with B3 + triplet memory.

The resulting y - B4_oof residual is therefore produced without using the outer row's label at any lower or triplet layer.

The final quartet memory is learned from these outer OOF B4 residuals.

### Hidden scoring

For hidden candidates:

- B2 and B3 use a base model fit on all visible rows.
- B4 adds triplet memory learned from cross-fitted visible residuals.
- B5 adds quartet memory learned from strict outer-OOF B4 residuals.

No hidden fitness value enters fitting or scoring.

## Higher-order rules

- Triplet minimum support: 2
- Triplet shrinkage: n/(n+2)
- Quartet minimum support: 2
- Quartet shrinkage: n/(n+4)
- Adaptive contribution is confidence-weighted and scaled by supported-subset coverage.
- No exact full-sequence memorization term.
- Candidate scope: 3 to 5 substitutions.

## Controls

The complete V8.1 cross-fit pipeline is rerun under:

1. GLOBAL_PERMUTED_B5 — global deterministic permutation of visible labels, seed 161.
2. COUNT_STRATIFIED_PERMUTED_B5 — deterministic permutation within mutation-count strata, seed 161.

## What would count as progress

This is not a pass/fail freeze. We are looking for a stable architecture candidate.

A breakpoint candidate should show:

- B5 remains above B3 in hidden ranking after cross-fitting.
- The gain is not confined to one tiny subgroup.
- Candidate selection changes meaningfully rather than only rescaling the same ranking.
- Both negative controls collapse.
- The same code can then be run unchanged on a second development landscape.

Only after replication do we create a canonical frozen architecture and a fresh sealed benchmark.
