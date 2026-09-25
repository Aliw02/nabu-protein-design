# Phase 2B — Candidate Assembly V1 Preregistration

## Goal

Move from selecting rows in a supplied candidate table to constructing valid mutation combinations from a known mutation vocabulary.

This stage tests **assembly only**.

It does not change the 2A acquisition controller and does not combine assembly with closed-loop acquisition. That combined experiment belongs to 2C.

## Required layers

### 1. ReferenceProtein

Every mutation must be checked against a reference protein.

For token `A10C`:

- position 10 must exist;
- the reference residue at position 10 must be `A`;
- target residue must be a standard amino acid;
- the mutation must change residue state.

A candidate cannot contain two substitutions at the same position.

### 2. MeasurementAggregator

Experimental replicate rows must be aggregated before entering frozen V8.3.

V1 frozen aggregation rule:

- canonicalize candidate identity;
- require finite labels;
- group exact canonical mutation sets;
- aggregate label by arithmetic mean;
- retain replicate count and sample standard deviation as metadata.

This layer is outside frozen V8.3.

### 3. CandidateAssembler V1

Inputs:

- fitted frozen V8.3 model;
- validated mutation vocabulary;
- measured mutation sets;
- target mutation order.

Search:

1. enumerate valid pairs from the vocabulary;
2. score pair candidates with frozen V8.3;
3. retain a bounded beam;
4. expand the beam by one valid mutation at a new position;
5. canonicalize, deduplicate and reject measured candidates;
6. score with frozen V8.3;
7. return top proposals by within-frontier V8.3 ordering.

No hidden fitness is used during assembly.

## Initial development POC

RhlA training pool is used only as a virtual post-freeze evaluator.

The assembler may see:

- known mutation vocabulary;
- bootstrap measured candidates and labels;
- reference identity information.

The assembler may not see:

- unmeasured RhlA labels;
- an unmeasured candidate ranking from the source table.

After proposals are frozen and hashed, the virtual evaluator may report which proposals happen to exist in the historical RhlA assay table and their known labels.

Unmatched proposals remain legitimate novel proposals with unknown virtual truth; they are not treated as failures.

## Initial V1 settings

- bootstrap: same deterministic 5% RhlA bootstrap;
- target assembly order: 3;
- pair beam width: 256;
- final proposal count: 64.

These are engineering search-budget settings, not biological optimum claims.

## Pass criteria

V1 passes engineering/scientific assembly POC when:

- all proposals satisfy reference constraints;
- no proposal duplicates a measured candidate;
- no duplicate proposal IDs exist;
- proposal list is deterministic;
- unmeasured truth perturbation cannot change pre-evaluation proposals;
- proposal IDs are frozen and hashed before virtual truth lookup;
- the search generates scoreable candidates where supported evidence exists.

Performance on historical RhlA labels is diagnostic only.
