# Phase 2C Freeze

## Status

**DEVELOPMENT COMPLETE — FROZEN WITH LIMITATION**

Frozen candidate for untouched Phase-2 evaluation:

`NABU_PHASE2C_GB1_GATED_V2`

## Frozen mechanism

The Phase-2C system uses:

- frozen Phase-1 V8.3 core;
- leak-safe closed-loop measurement;
- historical 50/50 acquisition before assembly maturity;
- CandidateAssembler V1;
- maturity gate based on repeated evidence;
- full CandidateAssembler beam probe;
- target-order triple assembly in the GB1 development POC;
- pre-reveal selection/proposal hashes.

Maturity rule:

- candidate outside measured set;
- candidate outside supplied reservoir;
- assayable identity;
- every main mutation support >= 2;
- every internal pair support >= 2;
- at least one complete batch of mature proposals.

The gate is evidence-driven, not percentage-driven.

## Engineering result

PASS:

- contract tests PASS;
- pinned GB1 source hash verified;
- hidden-truth perturbation invariant PASS in all four conditions;
- assembled measurements can re-enter the support model;
- deterministic shared bootstrap;
- proposal and acquisition lists frozen before reveal.

## Development scientific result

Preregistered best-fitness AUC gate was NOT passed.

Acquisition only:

- AUC = 0.691448;
- final best = 7.556869;
- final Top-1% hits = 3.

V2 acquisition + gated assembly:

- AUC = 0.691448;
- final best = 7.556869;
- final Top-1% hits = 10;
- assembly-source Top-1% hits = 7;
- assembly rounds = 2.

The combined system therefore increased elite-hit yield substantially without improving the preregistered best-so-far AUC.

## State-matched evidence

At 184 measurements, before truth reveal:

- assembly and acquisition batches were frozen from the exact same model state.

After reveal:

- assembly Top-1% hits = 4 / 8;
- acquisition Top-1% hits = 0 / 8;
- assembly mean fitness = 2.274018;
- acquisition mean fitness = 0.226603;
- assembly best = 5.430379;
- acquisition best = 1.615588.

This is strong mechanism evidence for maturity-gated assembly, but it does not retroactively change the preregistered V2 PASS rule.

## Freeze decision

Do not create V3 by further GB1 tuning.

Carry V2 forward unchanged as the Phase-2C frozen candidate, with the explicit limitation that its GB1 best-fitness AUC did not improve.

Final scientific judgment is deferred to the preregistered 2D blind benchmarks.

The AAV blind benchmark uses a separately preregistered Top-100 recovery primary metric and must not be altered based on GB1.
