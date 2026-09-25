# Phase 2C — Maturity-Gated Assembly + Acquisition Development Benchmark

## Goal

Run the first complete NABU closed-loop combinatorial-design benchmark in which measured evidence can mature into generation of candidates that were not supplied in the original acquisition reservoir.

This is a **real retrospective development evaluation**.

It is not a proof-of-concept and it is not the final blind Phase-2 verdict.

The final blind verdict remains quarantined in Phase 2D.

## Development landscapes

The complete protocol is run without per-landscape retuning on:

- GB1 Wu 2016 — reference four-site genotype `VDGV`;
- TrpB Johnston 2024 — reference four-site genotype `VFVS`;
- PhoQ — reference four-site genotype `AVST`.

These landscapes have already been used as development evidence and are permanently excluded from the final blind Phase-2 PASS.

## Reference-relative mutation coordinates

For every landscape, the four experimentally varied residues are mapped to benchmark-local positions 1..4.

A genotype is represented only by substitutions away from its biological four-site reference genotype.

No fitness label participates in mutation encoding or vocabulary construction.

## Assay universe versus supplied reservoir

For every landscape:

- **assay universe** = all measured genotype identities with finite historical fitness;
- **supplied reservoir** = deterministic identity-only subset of 2048 assay-universe genotypes;
- **bootstrap** = 64 deterministic identity-only genotypes from the reservoir;
- **final budget** = 256 total measurements;
- **batch size** = 8.

Acquisition-only methods may select only from the supplied reservoir.

Assembly may generate valid triple mutants outside the supplied reservoir.

Retrospective assay-universe membership may be used only as identity metadata so the virtual oracle can answer a proposal. Hidden fitness is not available to selection.

## Frozen acquisition references

- deterministic random;
- historical 50/50 support-deficit + V8.3 exploitation controller.

Relative-Evidence V2 is not added to this first 2C combination benchmark because 2A.4 did not establish it as superior to the historical 50/50 reference.

## CandidateAssembler V1 settings

Keep the frozen 2B mechanism:

- target order = 3;
- beam width = 256;
- proposal frontier retained for maturity filtering = 1024;
- reference-aware mutation validation;
- no measured duplicates.

No assembly parameter is changed between landscapes.

## Generic Evidence-Maturity Gate V1

At every round, CandidateAssembler V1 creates a pre-reveal proposal frontier.

A proposal is **mature** only if all of the following hold:

- valid against the biological four-site reference;
- outside the measured set;
- outside the supplied 2048-candidate reservoir;
- identity exists in the retrospective assay universe;
- every mutation has main-memory support >= 2;
- every internal mutation pair has pair-memory support >= 2.

The gate opens only when at least one complete experimental batch of 8 mature proposals exists.

Once opened, the gate is latched open.

Rationale:

- support >= 2 represents repeated measured evidence rather than an isolated observation;
- it is consistent with NABU's existing repeated-support higher-order memory convention;
- the rule is based on evidence state, not dataset identity or percentage measured.

No hidden fitness value is used by the gate.

## Conditions

Every condition uses the same assay universe, supplied reservoir, bootstrap, budget and batch cadence.

### 1. random_only

Deterministic random reservoir acquisition for the entire campaign.

### 2. acquisition_only

Historical 50/50 reservoir acquisition for the entire campaign.

### 3. random_then_gated_assembly

- deterministic random reservoir acquisition while the gate is closed;
- after the gate opens, experimental batches come only from mature assembled proposals.

This is the assembly control without intelligent acquisition warm-up.

### 4. acquisition_then_gated_assembly

- historical 50/50 reservoir acquisition while the gate is closed;
- after the gate opens, experimental batches come only from mature assembled proposals.

This is the primary maturity-gated combined system.

## State-matched first-gate ablation

At the first gate-open state of `acquisition_then_gated_assembly`:

- freeze the first mature assembly batch;
- independently freeze the historical-50/50 acquisition batch from the exact same pre-reveal model state;
- hash both lists before any truth lookup;
- evaluate the alternate acquisition batch only diagnostically;
- alternate labels never enter the main campaign.

This gives a direct same-state comparison of assembly versus acquisition.

## Anti-leakage requirements

Before every real campaign reveal:

- batch IDs are finalized;
- selection source is recorded;
- SHA256 of the batch is stored;
- maturity diagnostics are stored;
- no selected label has been returned by the oracle.

Selection code receives only:

- measured IDs and measured labels;
- assay-universe identities;
- supplied-reservoir identities;
- reference and mutation vocabulary;
- fitted V8.3 state.

It never receives unmeasured fitness values.

A hidden-truth perturbation test must show that changing unrevealed labels cannot alter the next pre-reveal decision from an identical measured state.

## Primary metrics

Per condition and landscape:

- post-bootstrap normalized discovery AUC;
- final best true fitness;
- final regret;
- cumulative Top-1% hits;
- measurements to first Top-1% discovery.

Assembly diagnostics:

- gate-open measurement;
- assembly-active rounds;
- number of assembled measurements;
- number of assembled measurements outside the original reservoir;
- assembled Top-1% hits;
- mature-proposal count by round;
- main/pair/triplet/quartet memory growth;
- router trajectory.

## Per-landscape scientific win criterion

The primary combined condition `acquisition_then_gated_assembly` is a development win over `acquisition_only` on a landscape only if:

1. post-bootstrap discovery AUC is strictly higher;
2. final regret is <= acquisition_only;
3. final Top-1% hit count is >= acquisition_only;
4. at least one assembly batch is actually measured.

## Overall 2C development verdict

- **PASS**: the combined system meets the per-landscape win criterion on at least two of the three landscapes and shows no clear joint degradation of both AUC and final regret on the remaining landscape.
- **PASS_WITH_LIMITATION**: at least one landscape meets the win criterion but the full suite is mixed.
- **NO_SCIENTIFIC_WIN**: no landscape meets the win criterion.
- **ENGINEERING_FAIL**: integrity, determinism, validity or leakage checks fail.

Do not alter the gate or assembly rule after seeing one landscape and then count later landscapes as preregistered evidence.

## Phase-2D boundary

No AAV/Random or FLIP2 IRED target labels may be loaded during 2C.

After 2C is complete, the entire system is frozen before the reserved external benchmarks are revealed.
