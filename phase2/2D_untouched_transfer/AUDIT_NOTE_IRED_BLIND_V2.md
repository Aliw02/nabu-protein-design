# Phase 2D IRED Blind V2 — Pre-Reveal Contract Failure Audit

Workflow run: `36165198162`  
Execution head: `eb31746803134fb73835b6093d0ea4201874c41d`

Status: **ENGINEERING CONTRACT FAILURE — NO BLIND DATA DOWNLOADED / NO TARGET REVEAL**

The workflow stopped at the synthetic Phase-2D contract tests. SciPy returned
Spearman `0.9999999999999999` for a perfect-ranking fixture while the test
asserted exact equality to `1.0`.

Fix:

- replace exact floating-point Spearman equality with `numpy.isclose`;
- no model logic changed;
- no metric definition changed;
- no source, split, scoreability, abstention policy, or target value was read
  to choose the fix.

The failed run remains part of the audit trail.
