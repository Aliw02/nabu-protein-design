# Phase 2C V2 Result

## Verdict

**Engineering PASS / Scientific development gate NOT PASSED.**

V2 changed only the maturity probe from 64 proposals to the full 256-proposal beam.

## Effect of the V2 change

The maturity gate opened earlier:

- V1 combined first gate: 208 measurements;
- V2 combined first gate: 184 measurements.

Assembly activation increased:

- V1 combined: 1 assembly round;
- V2 combined: 2 assembly rounds.

Final Top-1% discoveries:

- acquisition-only: 3;
- V2 combined: 10;
- assembly-source hits inside V2 combined: 7.

The full-beam gate therefore strengthened the assembly contribution.

## State-matched evidence

At 184 measurements, before any truth reveal:

- one assembly batch and one acquisition batch were both frozen from the same model state.

After freeze:

- assembly: 4 / 8 Top-1% hits;
- acquisition: 0 / 8 Top-1% hits;
- assembly best: 5.430379;
- acquisition best: 1.615588;
- assembly mean: 2.274018;
- acquisition mean: 0.226603.

This is strong evidence that the maturity gate identifies useful design states.

## Why the preregistered scientific gate still fails

The primary V1/V2 metric was best-fitness discovery AUC.

Acquisition-only had already found fitness 7.556869 before assembly activation.

The assembled high-fitness variants increased the number of elite discoveries but did not exceed that existing maximum, so best-so-far AUC remained exactly:

- acquisition-only: 0.691448;
- V2 combined: 0.691448.

The success criterion is not changed after observing this result.

## Decision

Do not create a V3 by repeatedly tuning the GB1 development landscape.

Before final 2D blind evaluation, perform a broader Phase-2C design review focused on objective definition:

- best-fitness optimization;
- elite-set discovery / hit yield;
- diversity of high-fitness discoveries.

The final blind benchmark's primary metric remains independently preregistered in FINAL_VALIDATION_PROTOCOL.md and is not changed by this GB1 result.
