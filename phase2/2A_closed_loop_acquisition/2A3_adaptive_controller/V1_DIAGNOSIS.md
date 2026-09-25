# 2A.3 V1 Diagnosis

## Outcome

Evidence-Adaptive V1 passed engineering checks but did not outperform the frozen historical 50/50 controller on the RhlA development landscape.

Primary micro_8 comparison:

- historical 50/50: Top-1% hits = 5, AUC = 0.976586;
- V1: Top-1% hits = 4, AUC = 0.969587.

Sensitivity micro_16:

- historical 50/50: Top-1% hits = 3, AUC = 0.970504;
- V1: Top-1% hits = 3, AUC = 0.953585.

Both V1 conditions still reached the global best.

## Controller trajectory diagnosis

Observed V1 pressure:

- initial exploration pressure: ~0.681;
- final exploration pressure: ~0.393;
- mean exploration pressure: ~0.49.

The B3/B4/B5 rank-disagreement gap was 0.0 throughout the summarized RhlA V1 trajectories.

Therefore V1 pressure was effectively driven by:

- structural support gap;
- scoreability gap.

The initial combined gap forced substantially more exploration than the historical 0.5 reference exactly when the historical controller was already discovering useful variants.

## V1 failure hypothesis

The absolute-gap formula confuses **remaining incompleteness** with **required exploration fraction**.

A support gap of ~0.59 does not imply that ~68% of acquisition pressure should be exploration.

This is a controller calibration problem, not evidence that support diagnostics are useless.

## V2 design constraint

Do not tune a replacement scalar from RhlA outcomes.

V2 must:

- preserve V1 as a recorded failure;
- use the bootstrap evidence state as its own internal reference;
- start neutral without a tuned constant;
- move exploration pressure according to relative evidence improvement or deterioration;
- keep the same V8.3 core and micro-batch cadence.
