# NABU V8.3 Canonical Phase-1 Core

**Status:** Frozen Phase-1 core architecture  
**Phase 1:** Closed  
**Core freeze:** `freeze/nabu-v8-3-dual-objective-router-validated`  
**Core commit:** `c8afdcd7698231d95a39ef13a3fe22b6e3f507c5`  
**Phase-1 closure freeze:** `freeze/nabu-phase1-final-cr9114-closed`  
**Closure commit:** `690494691c20edc3c9b7a6bba9a5156913173c58`

## Canonical rule

NABU V8.3 uses B3 raw-pair memory as the protected backbone and cross-fitted higher-order residual memory as a conditional correction.

The routing decision uses visible-data OOF behavior only.

1. If visible OOF B4 Spearman > B3 Spearman:
   - use global B5 higher-order scoring.
2. Otherwise inspect the B3-frozen top-20% elite region.
3. Allow B5 reranking inside that region only if:
   - OOF Top-50 mean true percentile improves, and
   - OOF Top-50 Top-1% hit count does not decrease.
4. Otherwise preserve B3 exactly.

Elite-only mode never allows higher-order scores to move a candidate across the B3-defined top-20% boundary.

## Final evidence by routing regime

### Global higher-order

Observed on:

- PHOT
- TrpB
- eqFP611
- CR9114-H1

CR9114-H1 final fresh gate:

- hidden 4/5 variants: 1,673
- B3 Spearman: 0.9271
- V8.3 Spearman: 0.9387
- strict gains: 4/5
- no regression: 5/5

### Elite-only higher-order

Observed on GB1.

Hidden result:

- B3 Spearman: 0.4018
- V8.3 Spearman: 0.4030
- Top-50 true score improved materially
- Top-50 Top-1% hits: 48 -> 50

### Protected B3

Observed when higher-order evidence is unsupported or harmful, including PhoQ and the low-order CreiLOV training regime.

PhoQ:

- B3 preserved at Spearman 0.5420
- harmful higher-order elite rerank was rejected

CreiLOV:

- B3-only extrapolation to 13,168 hidden 4/5-mutants
- Spearman 0.8968
- Top-1 regret 0.0333

## Active Acquisition is not part of the canonical core

The 50/50 exploration/exploitation controller is an optional outer policy.

Development freeze:

`freeze/nabu-v8-3-rhla-active-acquisition-dev-2x`

Evidence:

- RhlA retrospective Active 20% reached Passive 40% threshold.
- CR9114 fresh Active 10% reached Passive 40% threshold.
- CR9114 Active 20% failed the preregistered target.
- Active performance was non-monotonic across budgets.

Therefore Active Acquisition should be treated as a Phase-2 research component, not as part of the frozen V8.3 core.

## Preservation rule

Do not retune this Phase-1 core on:

- PHOT
- GB1
- TrpB
- PhoQ
- CreiLOV
- eqFP611
- RhlA
- CR9114-H1

Any Phase-2 modification should be implemented above or beside the frozen core and compared against it as a baseline.

## Phase-2 transition

The next stage changes the task from static ranking toward closed-loop combinatorial design.

See [PHASE2_PLAN.md](PHASE2_PLAN.md).
