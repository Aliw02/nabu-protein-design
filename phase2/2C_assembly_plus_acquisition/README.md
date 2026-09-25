# 2C — Assembly + Acquisition Development Benchmark

Phase 2C combines the validated Phase-2 acquisition loop with the reference-aware CandidateAssembler.

This stage is a **retrospective development benchmark**, not a POC and not the final blind evaluation.

## Development suite

The exact same protocol is applied to:

- GB1, reference genotype `VDGV`;
- TrpB, reference genotype `VFVS`;
- PhoQ, reference genotype `AVST`.

The controller and gate do not inspect the landscape name.

## Protocol

Per landscape:

- assay universe: all measured variants with finite fitness;
- supplied reservoir: 2048 identities selected without fitness;
- bootstrap: 64;
- final measurement budget: 256;
- batch size: 8;
- assembly target: triple mutants;
- beam width: 256;
- proposal frontier: 1024;
- maturity support threshold: 2.

## Maturity gate

Assembly does not activate at a fixed measurement percentage.

A proposal is mature only when:

- it is valid against the four-site reference;
- it is unmeasured;
- it is outside the supplied reservoir;
- it is assayable in the retrospective universe;
- every mutation has main support >= 2;
- every internal pair has pair support >= 2.

The gate opens when at least eight mature proposals exist and remains open thereafter.

## Conditions

- `random_only`
- `acquisition_only`
- `random_then_gated_assembly`
- `acquisition_then_gated_assembly`

The primary combined system is `acquisition_then_gated_assembly`.

At its first gate-open state, the benchmark freezes a same-state historical-50/50 acquisition batch and an assembly batch before either truth lookup, then reports their diagnostic outcomes separately.

## Integrity

Every selected batch is hashed before reveal.

The primary combined campaign is replayed independently on every landscape and its full selection transcript SHA256 must match.

## Final boundary

A 2C PASS remains development evidence.

AAV/Random and FLIP2 IRED remain quarantined until the entire 2C system is frozen for Phase 2D.
