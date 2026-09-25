# Phase 2C V1 Diagnosis

## Verdict

**Engineering PASS / Scientific development gate NOT PASSED.**

The preregistered primary condition did not strictly improve discovery AUC over acquisition-only:

- acquisition-only AUC: 0.691448;
- acquisition + gated assembly AUC: 0.691448.

Therefore V1 is not promoted as the frozen Phase-2C scientific design.

## Important positive mechanism evidence

The combined system increased final Top-1% discoveries:

- acquisition-only: 3;
- combined: 7;
- assembly-source Top-1% hits: 4.

At the first gate-open state (208 measurements), assembly and acquisition alternatives were frozen from the exact same pre-reveal model state.

State-matched result:

- assembly batch: 4 / 8 Top-1% hits;
- acquisition batch: 0 / 8 Top-1% hits;
- assembly best fitness: 5.430379;
- acquisition best fitness: 0.172249;
- assembly mean fitness: 2.401567;
- acquisition mean fitness: 0.035436.

This strongly supports that the maturity gate identified a state where assembled proposals were useful.

## Why AUC did not improve

The acquisition-only path had already found a high-fitness variant before assembly activated.

The preregistered AUC tracks **best-so-far fitness**, so additional high-fitness discoveries below that already-discovered maximum do not increase the curve.

This metric behavior is noted but the V1 success rule is not changed after seeing results.

## Gate behavior

With V1:

- historical acquisition first opened the gate at 208 measurements;
- the gate saw 9 mature candidates inside the top-64 proposal probe;
- one batch of 8 was measured;
- the next round had zero mature candidates inside the truncated top-64 probe;
- assembly therefore deactivated.

Random acquisition never opened the maturity gate.

This shows the gate reacts to evidence quality, but also exposes a design mismatch:

- CandidateAssembler beam width = 256;
- V1 maturity probe = only top 64.

Gate maturity was therefore judged on a truncated quarter of the actual search frontier.

## V2 single-change rationale

Do not alter:

- support >= 2 maturity rule;
- batch size;
- reservoir;
- bootstrap;
- total budget;
- acquisition rule;
- assembly beam width;
- primary metric;
- scientific success criteria.

Change only:

> maturity probe count: 64 -> 256

so the gate evaluates the full frozen CandidateAssembler beam frontier.

This is a structural consistency fix, not a label-fitted threshold change.

V1 remains permanently preserved as a non-pass.
