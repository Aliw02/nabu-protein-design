# Phase 2C — Maturity-Gated Assembly + Acquisition V1

## Goal

Test the first complete NABU closed-loop design system in which the model:

1. measures candidates from a supplied reservoir;
2. refits frozen V8.3;
3. checks whether evidence is mature enough for design;
4. when mature, assembles candidates that were not supplied in the reservoir;
5. freezes proposals before truth reveal;
6. assays them virtually;
7. learns from those measurements and repeats.

This is a development POC. It is not the final blind Phase-2 verdict.

## Development landscape

GB1 Wu 2016 four-site landscape.

Biological reference:

- Protein G B1 domain, 56 aa;
- mutable residues: V39, D40, G41, V54;
- WT four-site genotype: VDGV.

GB1 is already a development dataset and therefore is permanently excluded from the final blind Phase-2 PASS.

## Known candidate reservoir versus design universe

To test whether assembly adds something beyond a supplied candidate table:

- the **assay universe** is the available GB1 landscape identities;
- the **supplied candidate reservoir** is a deterministic identity-only 2048-candidate subset;
- every condition begins from the same 64-candidate identity-only bootstrap inside that reservoir;
- acquisition-only methods may select only from the supplied reservoir;
- assembly may propose valid GB1 triple mutants outside the supplied reservoir;
- proposal eligibility may use identity availability but never hidden fitness.

Final measurement budget:

- 256 total measurements.

Batch size:

- 8.

## Frozen acquisition references

- deterministic random;
- historical 50/50.

Relative-Evidence V2 is not introduced as another moving part in the first 2C experiment.

## Generic maturity gate

Do **not** open assembly at a fixed measurement percentage.

At every round, run the frozen CandidateAssembler V1 proposal probe.

A proposed triple is **mature** only when:

- it is valid against the GB1 reference;
- it is outside the measured set;
- it is outside the supplied candidate reservoir;
- it exists in the identity-only assay universe;
- every mutation has main-memory support >= 2;
- every internal mutation pair has pair-memory support >= 2.

Assembly opens only if the current model can supply at least one complete experimental batch of 8 mature proposals.

Rationale:

- support >= 2 means repeated measured evidence rather than one isolated observation;
- this matches the repeated-support convention already used by NABU higher-order memories;
- the gate depends on evidence structure, not dataset name or percent measured.

No target label is used by the gate.

## Conditions

All use the same bootstrap, reservoir, oracle and total budget.

1. **random_only**
   - random reservoir acquisition for the full campaign.

2. **acquisition_only**
   - historical 50/50 reservoir acquisition for the full campaign.

3. **random_then_gated_assembly**
   - random acquisition while gate is closed;
   - once gate opens, use mature assembled proposals.

4. **acquisition_then_gated_assembly**
   - historical 50/50 acquisition while gate is closed;
   - once gate opens, use mature assembled proposals.

## State-matched assembly ablation

At the first gate-open state of the combined campaign:

- freeze the assembly batch;
- freeze the historical-50/50 acquisition batch from the exact same pre-reveal model state;
- only then look up both batches' virtual truth for diagnostic comparison;
- the alternate acquisition labels must not affect the main campaign.

This isolates immediate assembly quality from state differences.

## Anti-leakage requirements

- hidden labels are accessible only through the oracle;
- every batch ID list is hashed before reveal;
- gate diagnostics are recorded before reveal;
- assembled proposals are frozen before reveal;
- changing all currently unrevealed labels while preserving measured labels must not change the next decision;
- no blind AAV or IRED benchmark labels may be loaded.

## Primary metric

Post-bootstrap normalized discovery AUC over the full GB1 assay universe.

## Scientific 2C V1 gate

Engineering PASS requires:

- deterministic replay;
- hidden-truth perturbation invariance;
- pre-reveal hashes;
- no invalid assembled candidates;
- no duplicate measurements;
- if assembly activates, every assembled candidate is outside the supplied reservoir and passes the maturity rule.

Scientific development PASS requires the primary combined condition
`acquisition_then_gated_assembly` to satisfy all of:

1. discovery AUC > acquisition_only;
2. final regret <= acquisition_only;
3. final Top-1% hit count >= acquisition_only;
4. at least one post-gate assembly batch is actually measured.

If these are not met, preserve the result as a 2C development failure/limitation rather than tuning against the same run.

## Boundary

A 2C development PASS does not equal final Phase-2 PASS.

After 2C is frozen, final evidence comes only from the reserved 2D blind benchmarks defined in `phase2/FINAL_VALIDATION_PROTOCOL.md`.
