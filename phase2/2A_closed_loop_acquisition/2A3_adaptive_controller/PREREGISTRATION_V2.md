# Phase 2A.3 — Relative-Evidence Controller V2 Preregistration

## Motivation

V1 used absolute evidence gaps directly and over-explored early on RhlA.

V2 changes the **pressure calibration rule**, not the evidence sources or frozen core.

## Evidence gap

Use the same state gaps:

- support gap;
- scoreability gap;
- hierarchy rank-disagreement gap.

Combine them as:

```
current_gap =
    1 - (
        (1 - support_gap)
        * (1 - scoreability_gap)
        * (1 - disagreement_gap)
    )
```

## Self-reference calibration

At the first post-bootstrap acquisition round, store:

`reference_gap = current_gap`

Then compute:

```
exploration_pressure =
    current_gap / (current_gap + reference_gap)
```

Properties:

- at bootstrap: pressure = 0.5 exactly;
- if evidence gaps shrink: pressure falls below 0.5;
- if evidence gaps worsen: pressure rises above 0.5;
- no fitted scalar or dataset-specific threshold is introduced.

If both gaps are zero, pressure is zero.

## Candidate exploration evidence

Use:

- structural exploration percentile;
- hierarchy rank disagreement.

Candidate exploration evidence is their maximum.

Unscoreable candidates retain zero exploitation rank and can still be selected through structural exploration.

## Comparison

Primary:

- historical 50/50 + micro_8;
- relative_evidence_v2 + micro_8.

Sensitivity:

- historical 50/50 + micro_16;
- relative_evidence_v2 + micro_16.

Same 5% bootstrap and 20% total budget.

## Interpretation boundary

V2 is developed after observing V1 failure on RhlA.

Therefore RhlA V2 results are development evidence only.

Generalization must be tested in 2A.4 on multiple landscapes without per-landscape retuning.
