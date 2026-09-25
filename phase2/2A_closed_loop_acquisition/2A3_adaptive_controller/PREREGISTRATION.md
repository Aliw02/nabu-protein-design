# Phase 2A.3 — Evidence-Adaptive Controller V1 Preregistration

## Question

Can a controller that changes exploration pressure from the current evidence state improve or preserve discovery relative to the frozen historical 50/50 controller?

## Single changed mechanism

The frozen V8.3 core is unchanged.

The campaign simulator is unchanged.

The primary cadence is frozen from 2A.2:

- primary: micro_8;
- sensitivity: micro_16.

The only new mechanism is the acquisition controller.

## Reference baseline

`historical_50_50`

The baseline definition from 2A.1 remains unchanged.

## Evidence-Adaptive V1

V1 uses only information available before label reveal.

For each round, compute three state-level knowledge gaps:

1. **support gap**
   - mean historical structural support deficit across the unmeasured pool;
   - main and pair support only.

2. **scoreability gap**
   - fraction of the unmeasured pool not currently scoreable by V8.3.

3. **rank-disagreement gap**
   - mean within-candidate disagreement among B3/B4/B5 percentile ranks;
   - zero when the hierarchy agrees.

These are combined without tuned scalar weights:

```
knowledge_completeness =
    (1 - support_gap)
    * (1 - scoreability_gap)
    * (1 - disagreement_gap)

exploration_pressure =
    1 - knowledge_completeness
```

Therefore:

- sparse support increases exploration;
- poor scoreability increases exploration;
- hierarchy disagreement increases exploration;
- when all three improve, exploitation pressure rises automatically.

No dataset name appears in the rule.

## Candidate evidence

Each candidate receives:

- structural exploration percentile;
- B3/B4/B5 rank disagreement;
- unscoreable indicator;
- current V8.3 exploitation percentile.

Candidate exploration evidence is the maximum of:

- structural exploration percentile;
- hierarchy disagreement;
- unscoreable indicator.

Final acquisition score:

```
exploration_pressure * exploration_evidence
+ (1 - exploration_pressure) * exploitation_rank
```

The pressure is recomputed after every reveal/refit.

## Important boundary

This is **not** calibrated epistemic uncertainty.

Support and cross-hierarchy disagreement are evidence diagnostics only.

## Experiment matrix

Run four conditions to 20% budget from the same 5% bootstrap:

1. historical_50_50 + micro_8
2. evidence_adaptive_v1 + micro_8
3. historical_50_50 + micro_16
4. evidence_adaptive_v1 + micro_16

## Primary metrics

- normalized post-bootstrap discovery AUC;
- final best true fitness;
- final regret;
- cumulative Top-1% hits;
- first Top-1% measurement.

Diagnostics:

- scoreable coverage;
- B3 OOF Spearman;
- router trajectory;
- exploration-pressure trajectory.

## Interpretation rule

RhlA is a development landscape.

This experiment may determine whether V1 is worth carrying into 2A.4, but it cannot establish general superiority.

If V1 fails, preserve the failure and diagnose which evidence term caused the behavior before trying a V2.
