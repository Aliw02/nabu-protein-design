# Phase 2C Freeze — V3 Dynamic Maturity Gate

## Status

**PASS — development benchmark complete.**

Authoritative Phase-2C version:

`NABU_PHASE2C_MULTILANDSCAPE_V3`

CI validation:

- workflow run: `36160574928`;
- validated head: `316910cfc96a4bb98027b01c6f7cbba91646edc2`;
- contract tests: **17 passed**;
- benchmark job: PASS;
- artifact: `10876390166`;
- artifact digest: `sha256:d12852968cd865775e34abd3474f04e4226e98440b70b6093e8f48763b5c152d`.

Earlier GB1-only V1/V2 attempts remain in repository history and are not deleted.

The multilandscape V1 latch failure is preserved in `V1_FAILURE.md`.

## Frozen V3 architecture

Per round:

1. fit frozen V8.3 on measured data;
2. evaluate current mature triple-mutant assembly frontier;
3. if >= one full batch of 8 mature proposals exists, select an assembly batch;
4. otherwise use the condition's frozen reservoir acquisition policy;
5. freeze candidate IDs and SHA256 before label reveal;
6. virtual assay selected candidates only;
7. update measured state and refit.

Mature proposal rule:

- reference-valid;
- unmeasured;
- outside the original 2048-candidate reservoir;
- assayable in the retrospective universe;
- every mutation main support >= 2;
- every internal pair support >= 2.

Frozen search settings:

- target order = 3;
- beam width = 256;
- proposal frontier = 1024;
- batch size = 8.

The gate is dynamic per round. It is not a fixed measurement-percentage threshold.

## Development benchmark

Identical protocol was run on:

- GB1 / reference VDGV;
- TrpB / reference VFVS;
- PhoQ / reference AVST.

Each used:

- deterministic identity-only reservoir = 2048;
- shared bootstrap = 64;
- final budget = 256.

Conditions:

- random_only;
- acquisition_only;
- random_then_gated_assembly;
- acquisition_then_gated_assembly.

## Primary combined result

The preregistered primary condition was:

`acquisition_then_gated_assembly`

Compared with `acquisition_only`:

| Landscape | Δ discovery AUC | Δ regret | Δ Top-1% hits | Combined win |
| --- | ---: | ---: | ---: | --- |
| GB1 | +0.071326 | -1.999854 | +6 | YES |
| TrpB | +0.038913 | -0.065887 | +10 | YES |
| PhoQ | +0.026854 | -4.265920 | +4 | YES |

No landscape showed the preregistered joint degradation of both lower AUC and higher final regret.

Therefore the frozen overall 2C development verdict is:

> **PASS**

## Absolute combined outcomes

GB1 combined:

- discovery AUC = 0.417691;
- final best = 5.063573;
- final Top-1% hits = 9;
- first assembly-ready state = measurement 112;
- assembly measurements = 56;
- assembled Top-1% hits = 6.

TrpB combined:

- discovery AUC = 0.568869;
- final best = 0.596605;
- final Top-1% hits = 12;
- first assembly-ready state = measurement 120;
- assembly measurements = 48;
- assembled Top-1% hits = 10.

PhoQ combined:

- discovery AUC = 0.142330;
- final best = 27.112305;
- final Top-1% hits = 7;
- first assembly-ready state = measurement 80;
- assembly measurements = 48;
- assembled Top-1% hits = 4.

Every assembly-selected candidate in the primary combined runs was outside the original supplied reservoir.

## Same-state first-ready diagnostic

At the first ready model state, assembly and historical-50/50 acquisition batches were both frozen before truth lookup.

GB1 at measurement 112:

- assembly batch mean fitness = 0.664056;
- assembly best = 2.454478;
- assembly Top-1% hits = 1;
- acquisition batch mean fitness = 0.071697;
- acquisition best = 0.541482;
- acquisition Top-1% hits = 0.

TrpB at measurement 120:

- assembly batch mean fitness = 0.250883;
- assembly best = 0.596605;
- assembly Top-1% hits = 2;
- acquisition batch mean fitness = 0.016748;
- acquisition best = 0.030198;
- acquisition Top-1% hits = 0.

PhoQ at measurement 80:

- assembly batch mean fitness = 1.382275;
- assembly best = 6.902132;
- assembly Top-1% hits = 0;
- acquisition batch mean fitness = 0.002163;
- acquisition best = 0.014637;
- acquisition Top-1% hits = 0.

This diagnostic supports the immediate usefulness of mature assembly at the exact same pre-reveal model state on all three development landscapes.

## Gate behavior

Maturity was intermittent rather than permanently stable.

Primary combined condition:

- GB1: 7 assembly-ready rounds; 11 post-first-readiness fallback rounds; 12 ready/not-ready transitions.
- TrpB: 6 assembly-ready rounds; 11 post-first-readiness fallback rounds; 12 transitions.
- PhoQ: 6 assembly-ready rounds; 16 post-first-readiness fallback rounds; 12 transitions.

This validates why the V1 permanently latched gate was structurally incorrect and why V3 uses per-round readiness.

## Integrity checks

PASS:

- deterministic full-campaign replay on GB1, TrpB and PhoQ;
- hidden/unrevealed truth perturbation invariant for the first pre-reveal decision on all three landscapes;
- source SHA256 verification;
- pre-reveal selection hashes;
- no duplicate measurements;
- reference-aware mutation encoding;
- reserved Phase-2D blind datasets were not loaded.

## Scientific boundary

This is a strong development PASS, not the final Phase-2 scientific PASS.

GB1, TrpB and PhoQ have already influenced development.

Do not retune V3 on these results before the reserved external test.

## Next gate

Phase 2C is frozen.

The next stage is:

**2D — Untouched External Benchmark**

Only now may the reserved AAV/Random and FLIP2 IRED benchmark labels be revealed under `phase2/FINAL_VALIDATION_PROTOCOL.md`.
