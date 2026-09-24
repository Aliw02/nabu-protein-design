# NABU Phase 2 — Closed-Loop Combinatorial Design

## Goal

Phase 2 changes the task.

Phase 1 asked:

> Given a candidate library, can NABU learn the landscape and rank unseen variants?

Phase 2 asks:

> Given a mutation vocabulary and a limited experimental budget, can NABU actively assemble, choose, and iteratively discover high-fitness variants?

This is the first transition from static landscape inference toward an actual combinatorial protein design loop.

## Important boundary

Phase 2 is **not yet full de novo protein generation**.

The initial design space remains a known protein with a defined mutation vocabulary or mutable positions.

NABU will be allowed to propose new combinations inside that space, including combinations not present in the measured training set.

Full sequence generation from scratch, structural generation, or AlphaFold-like structure inference are later problems.

## Phase 2A — closed-loop acquisition

Start from a small measured seed set.

For each round:

1. fit the frozen V8.3 core,
2. evaluate current interaction memory and uncertainty/support structure,
3. select a new batch to measure,
4. reveal only that batch's virtual-assay truth,
5. update memory,
6. repeat.

Primary comparison:

- random acquisition,
- greedy exploitation,
- exploration-only,
- current 50/50 controller,
- new Phase-2 controller candidates,
- frozen V8.3 static baseline.

Primary metric is no longer only Spearman.

We care about:

- best true fitness discovered after N measurements,
- experiments required to first reach Top-1%,
- cumulative Top-1% discoveries,
- regret over acquisition rounds,
- diversity of discovered solutions,
- measurement efficiency.

## Phase 2B — combinatorial candidate assembly

Instead of choosing only from a supplied candidate table, allow NABU to construct mutation combinations from the known mutation vocabulary.

Example:

```text
known mutation vocabulary
        ↓
main / pair / triplet / quartet memories
        ↓
assemble candidate combinations
        ↓
constraint and support filters
        ↓
rank proposed designs
```

The first implementation should remain conservative:

- no unknown amino-acid vocabulary,
- no full-sequence generation,
- no hardcoded biological wordlists,
- no structural model requirement,
- avoid combinatorial explosion using adaptive search rather than fixed hop limits.

## Phase 2C — iterative design campaign

Combine candidate assembly with active measurement.

The loop becomes:

```text
measured seed set
    ↓
fit NABU
    ↓
assemble novel combinations
    ↓
choose experimental batch
    ↓
virtual assay
    ↓
update memory
    ↓
repeat
```

Use a fixed total measurement budget so different policies can be compared fairly.

## Phase 2D — untouched transfer

After controller/design development is complete, freeze it.

Then run one new untouched landscape as the Phase-2 validation.

No tuning after reveal.

## Core preservation

The Phase-1 V8.3 core remains the baseline:

`freeze/nabu-v8-3-dual-objective-router-validated`

Any Phase-2 improvement should be implemented as an outer design/search/acquisition layer or as an explicitly versioned new architecture.

The Phase-1 result must remain reproducible and untouched.

## First implementation target

The first Phase-2 experiment should be a **virtual wet-lab campaign** on a complete combinatorial landscape:

- begin with approximately 5% measured variants,
- use multiple acquisition rounds,
- allow generation/assembly of unseen combinations,
- reveal only selected candidate truth each round,
- compare discovery efficiency against random and static ranking.

Success means NABU discovers high-fitness variants with fewer measurements, not merely that it obtains a higher global correlation.
