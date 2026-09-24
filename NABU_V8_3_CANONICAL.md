# NABU V8.3 Canonical Working Architecture

**Status:** Current canonical development architecture  
**Date:** 2026-09-24  
**Base:** V8.1 cross-fitted higher-order hierarchy  
**Supersedes as current selector:** V8.2 single-objective OOF gate and the first V8.3 boundary-only router

## Architecture

NABU V8.3 uses B3 raw-pair memory as the protected backbone and cross-fitted higher-order residual memory as a conditional correction.

The router is determined only from visible-data OOF behavior:

1. If visible OOF B4 Spearman is higher than B3 Spearman, use global B5 higher-order scoring.
2. Otherwise inspect elite OOF behavior inside the B3-frozen top-20% region.
3. Allow rank-preserving B5 reranking inside that region only if:
   - OOF Top-50 mean true percentile improves, and
   - OOF Top-50 Top-1% hit count does not decrease.
4. Otherwise protect B3 and apply no higher-order reranking.

For elite-only mode, B3 freezes region membership. Higher-order is allowed to reorder candidates only inside the frozen region, so no cross-boundary score distortion is possible.

## Current evidence

### PHOT

Mode: GLOBAL_HIGHER_ORDER

- B3 Spearman: 0.899511
- V8.3 / B5 Spearman: 0.911189
- Top-50 Top-1% hits: 10 -> 11
- Top-50 mean true score: 1.246918 -> 1.264681

The dual-objective extension does not change the positive-global-OOF branch, so PHOT retains the established V8.3 result.

### GB1

Mode: RANK_PRESERVING_B3_TOP20_B5_RERANK

Visible OOF:
- global B4-B3 Spearman delta: -0.002357
- elite Top-50 percentile delta: +0.000319
- elite Top-1% hit delta: 0

Hidden:
- B3 Spearman: 0.401773
- V8.3 Spearman: 0.403033
- Top-50 Top-1% hits: 48 -> 50
- Top-50 mean true score: 4.320282 -> 4.520412
- six promoted Top-50 variants mean true score: 4.692752
- six displaced variants mean true score: 3.025001

### TrpB Johnston 2024

Mode: GLOBAL_HIGHER_ORDER

Visible OOF:
- global B4-B3 Spearman delta: +0.007429
- elite Top-50 percentile delta: +0.001382
- elite Top-1% hit delta: +3

Hidden:
- B3 Spearman: 0.290170
- V8.3 Spearman: 0.299617
- Top-50 Top-1% hits: 45 -> 47
- Top-50 mean true score: 0.584520 -> 0.615493
- normalized Top-1 regret: 0.318988 -> 0.034064
- seven promoted Top-50 variants mean true score exceeded displaced variants by +0.221234

### PhoQ

Mode: B3_PROTECTED_NO_HIGHER_ORDER

Visible OOF:
- global B4-B3 Spearman delta: -0.002613
- elite Top-50 percentile delta: +0.000036
- elite Top-1% hit delta: -4

The elite safety gate therefore rejected higher-order reranking.

Hidden:
- B3 Spearman preserved: 0.541983
- Top-50 Top-1% hits preserved: 31
- Top-50 mean true score preserved: 17.373217

Without the elite safety gate, higher-order would have reduced hidden Top-50 hits from 31 to 29 and mean Top-50 true score from 17.373217 to 16.304157.

## Current conclusion

Across the four development landscapes currently available:

- PHOT: higher-order is useful globally.
- TrpB: higher-order is useful globally.
- GB1: higher-order is useful for elite reranking but harmful globally.
- PhoQ: higher-order should be suppressed.

The useful result is therefore not a universal higher-order model. It is a visible-OOF-routed hierarchy that can choose between global higher-order, elite-only higher-order, and protected pairwise scoring.

This is the current NABU architecture to preserve. Future validation should test this frozen rule without changing its routing thresholds or decision logic.
