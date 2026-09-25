# Phase 2A.2 Development Freeze

## Status

**PASS — batch staleness supported on the RhlA development landscape.**

## Experiment invariant

The acquisition rule was fixed:

`historical_50_50`

Only reveal/refit cadence changed.

All conditions shared:

- the same 5% identity-only bootstrap;
- the same candidate pool;
- the same virtual assay oracle;
- the same frozen V8.3 core;
- the same 20% final measurement budget;
- the same discovery metrics.

## Observed comparison

| Cadence | Refits | Max batch | Best fitness | Regret | Top-1% hits | First Top-1% | Discovery AUC |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| historical_jumps | 2 | 161 | -0.010777 | 0.717198 | 2 | 322 | 0.922025 |
| micro_32 | 8 | 32 | 0.706421 | 0.000000 | 3 | 273 | 0.929843 |
| micro_16 | 16 | 16 | 0.706421 | 0.000000 | 3 | 113 | 0.970504 |
| micro_8 | 31 | 8 | 0.706421 | 0.000000 | 5 | 105 | 0.976586 |
| micro_4 | 61 | 4 | 0.706421 | 0.000000 | 3 | 193 | 0.952646 |

## Scientific interpretation

Supported on this development landscape:

- large stale batches can severely degrade acquisition;
- micro-batching materially improves discovery under the same 50/50 controller;
- smaller batches are not monotonically better;
- micro_8 produced the strongest observed RhlA result;
- micro_16 was also strong and is retained as a sensitivity control.

Not established:

- a universal optimal batch size;
- superiority of micro_8 on unseen protein landscapes;
- any new claim about B4/B5, because final router mode remained B3-protected in the summarized conditions.

## Freeze decision

The 2A.2 experiment definition and interpretation are frozen for the next stage.

For 2A.3:

- use micro_8 as the primary development cadence;
- repeat key controller comparisons with micro_16 as a sensitivity check;
- change controller logic only;
- preserve the frozen historical 50/50 controller as the baseline.

Do not reopen 2A.2 by tuning batch size against later 2A.3 outcomes.
