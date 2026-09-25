# Phase 2A.3 Development Freeze

## Status

**DEVELOPMENT COMPLETE — carry V2 forward as an adaptive candidate, with historical 50/50 retained as the reference.**

## V1 outcome

Evidence-Adaptive V1 was engineered successfully but underperformed the historical 50/50 controller on RhlA.

V1 is preserved as a failure and is not rewritten.

Diagnosis:

- absolute evidence gaps produced excessive early exploration;
- initial exploration pressure was ~0.681;
- hierarchy disagreement was 0.0 throughout the summarized RhlA trajectory;
- support and scoreability gaps dominated the controller.

See `V1_DIAGNOSIS.md`.

## V2 change

V2 uses relative evidence calibration.

At bootstrap:

- exploration pressure = 0.5 exactly.

As evidence gaps shrink:

- exploration pressure falls below 0.5.

No fitted scalar or dataset-name threshold is used.

Observed RhlA trajectory:

- micro_8: pressure 0.500 → 0.367;
- micro_16: pressure 0.500 → 0.364.

The controller therefore adapts in the intended direction.

## RhlA development results

| Condition | Best | Regret | Top-1% hits | First Top-1% | AUC |
| --- | ---: | ---: | ---: | ---: | ---: |
| historical 50/50, micro_8 | 0.706421 | 0.000000 | 5 | 105 | 0.976586 |
| relative V2, micro_8 | 0.706421 | 0.000000 | 4 | 105 | 0.971920 |
| historical 50/50, micro_16 | 0.706421 | 0.000000 | 3 | 113 | 0.970504 |
| relative V2, micro_16 | 0.706421 | 0.000000 | 4 | 113 | 0.968171 |

## Interpretation

V2 is not a universal winner on RhlA.

It improved materially over V1 and behaved as designed, but:

- micro_8 historical 50/50 retained higher AUC and one more Top-1% hit;
- micro_16 V2 found one additional Top-1% hit but had slightly lower AUC.

This is sufficient to define a stable adaptive candidate for generalization testing, but not sufficient to claim superiority.

## Freeze decision for 2A.4

Carry forward exactly these controllers:

1. `historical_50_50` — fixed reference;
2. `relative_evidence_v2` — adaptive candidate.

Carry forward cadence:

- primary: micro_8;
- sensitivity: micro_16 where computationally practical.

Do not create V3 from further RhlA tuning before multi-landscape evaluation.

2A.4 must determine whether V2 generalizes across multiple landscapes without per-landscape retuning.
