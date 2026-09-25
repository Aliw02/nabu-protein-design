# 2C — Assembly + Acquisition

This stage combines the two Phase-2 mechanisms for the first time.

Development POC:

- GB1 Wu 2016;
- true GB1 reference sequence;
- mutable sites V39, D40, G41 and V54;
- supplied identity-only reservoir of 2048 variants;
- 64-measurement shared bootstrap;
- 256-measurement final budget;
- batch size 8.

## Maturity gate

Assembly is not enabled by a fixed measurement percentage.

The gate opens only when CandidateAssembler V1 can provide a complete batch of novel-to-reservoir, assayable triple-mutant proposals for which:

- every main mutation has support >= 2;
- every internal pair has support >= 2.

## Controls

- random only;
- historical 50/50 acquisition only;
- random then gated assembly;
- historical 50/50 acquisition then gated assembly.

At the first gate-open state, a state-matched assembly-vs-acquisition batch ablation is also frozen before virtual truth lookup.

## Final boundary

GB1 is development-only evidence.

The reserved blind AAV and IRED benchmarks remain quarantined until 2C is frozen.
