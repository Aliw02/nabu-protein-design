# Phase 2C V1 Failure Record

## Status

**ENGINEERING FAIL — preserved.**

Workflow:

- run: 36159740144
- head: 6f24a7ab39fd446aa22189aa9a915ff8d9b74a3a
- artifact: 10874338883
- artifact digest: sha256:7b6e0315474e5f5ada583478e0b17c7e613732e74c605c72f27b26683cc543ca

## What passed before failure

- contract tests;
- pinned source downloads;
- source SHA256 verification;
- GB1 random_only completed;
- GB1 acquisition_only completed.

Observed control results before failure:

- GB1 random_only: best 3.852525, Top-1% hits 4, discovery AUC 0.386663;
- GB1 acquisition_only: best 3.063719, Top-1% hits 3, discovery AUC 0.346365.

These partial values are retained but are not a completed 2C comparison.

## Failure

During `GB1/random_then_gated_assembly`:

- the V1 maturity gate opened after finding at least one complete mature batch;
- after an assembly batch was measured and the model refit, the next assembly frontier contained only 4 mature proposals;
- V1 had latched the gate permanently open;
- the campaign therefore could not supply the required batch of 8 and stopped.

Runtime error:

`latched assembly gate cannot supply batch of 8; mature=4`

## Diagnosis

The scientific support rule did not fail.

The failure was caused by the **persistence assumption**:

> once assembly is mature once, it will remain able to supply a complete batch forever.

That assumption is not guaranteed because:

- selected mature proposals are removed after measurement;
- the ranked beam/frontier changes after refit;
- support can increase while the currently retained mature proposal count still falls.

## Allowed V2 change

V2 may change only the gate persistence rule.

Keep frozen:

- support threshold = 2;
- target order = 3;
- beam width = 256;
- proposal frontier = 1024;
- bootstrap = 64;
- total budget = 256;
- batch size = 8;
- all conditions;
- all metrics;
- all scientific win criteria.

V1 results and failure remain preserved.
