# Phase 2D IRED — First Blind Reveal Audit

Date: 2026-09-26
Status: PRESERVED_FIRST_BLIND_REVEAL

## Dataset

FLIP2 IRED two-to-many

Source SHA256:
aa45a2f85fb1af87b6b0e86397b3f8f292a061dc574be2536673b5bd490b0e74

Fit rows: 3746
Validation rows: 662
Official test rows: 4178

No Phase-2D scientific tuning was performed after target reveal.

## Pre-reveal structural result

Strict frozen V8.3 scoreable test rows:
131 / 4178 = 3.1354715%

Unscoreable:
4047 / 4178

Rows with all main effects supported:
3648 / 4178

## Actual policy used in this first reveal

IRED_ABSTENTION_V2

Unscoreable rows were assigned an abstention floor below all scoreable
frozen-router predictions.

B2 fallback used: false

This first blind result is preserved exactly as produced and must not be
overwritten by any later audited rerun.

## Frozen router result

Router mode:
B3_PROTECTED_NO_HIGHER_ORDER

Visible OOF Spearman:
0.32860921025867007

Triplet memory entries:
0

Quartet memory entries:
0

## Blind full-test result

Spearman:
-0.021770573158261506

NDCG:
0.9426983698477694

Top-1% K:
42

Top-1% hits:
0

Top-1% recall:
0.0

Top-1% enrichment:
0.0

Normalized regret:
0.23152235169646065

Best true value among predicted Top-1%:
1.1889779702946226

## Strict-scoreable diagnostic only

Count:
131

Spearman:
0.11938824534244381

NDCG:
0.8771807667453829

This subset diagnostic is not the primary Phase-2D result.

## Audit rule

Any subsequent change to full-test handling must:

1. preserve this directory unchanged;
2. document the reason separately;
3. use a new version;
4. preserve the new result separately;
5. never retroactively replace this first blind reveal.
