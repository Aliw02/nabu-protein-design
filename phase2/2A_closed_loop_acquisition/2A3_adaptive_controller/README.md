# 2A.3 — Evidence-Adaptive Controller

This stage changes the acquisition controller while keeping the frozen Phase-1 V8.3 core and the 2A.2 cadence decision fixed.

## V1 principle

Exploration pressure is derived from the current evidence state:

- structural support gap;
- scoreability gap;
- B3/B4/B5 rank disagreement.

No dataset-name rule is used.

No fixed 50/50 controller weight is tuned.

## Development protocol

Primary cadence:

- micro_8.

Sensitivity cadence:

- micro_16.

Reference:

- historical 50/50.

Outputs:

- per-condition discovery curves;
- adaptive controller pressure trajectories;
- controller comparison;
- controller manifest.

Read `PREREGISTRATION.md` before interpreting results.
