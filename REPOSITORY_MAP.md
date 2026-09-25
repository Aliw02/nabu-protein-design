# NABU Repository Map

This document separates the repository by role so the scientific audit trail and the reusable software package are not confused.

## 1. Canonical reusable runtime

```text
src/nabu_protein/
    core.py           Legacy pairwise package API kept for compatibility
    higher_order.py   Frozen Phase-1 B2/B3/B4/B5 hierarchy
    router.py         Frozen V8.3 dual-objective OOF router
    v83.py            Reusable facade over the frozen V8.3 core
    cli.py            Existing legacy CLI; not yet migrated to the V8.3 facade
```

The canonical Phase-1 scientific source is:

```text
freeze/nabu-v8-3-dual-objective-router-validated
c8afdcd7698231d95a39ef13a3fe22b6e3f507c5
```

The package hardening work must reproduce that implementation rather than reinterpret it.

## 2. Historical scientific experiment code

```text
experiments/
    v7_higher_order/
    v8_adaptive_higher_order_dev/
    v8_1_crossfit_higher_order_dev/
    v8_1_gb1_replication_dev/
    v8_2_objective_tournament/
    v8_3_adaptive_router/
    v8_3_multilandscape_validation/
    v8_3_creilov_sealed/
    v8_3_eqfp611_sealed/
    v8_3_rhla_sample_efficiency/
    v8_3_rhla_active_acquisition/
    v8_3_cr9114_final_phase1/
```

These directories are historical scientific records and parity references. Do not rewrite them as part of package cleanup.

## 3. Sealed inputs

Current sealed-input roots include:

```text
creilov_sealed_input/
eqfp611_sealed_input/
rhla_sample_efficiency_sealed_input/
cr9114_h1_final_gate_sealed_input/
```

Their manifests and hashes are part of the audit trail. They are not package source code.

## 4. Frozen/result artifacts

Top-level `nabu_v8_*` result directories are historical outputs from preregistered or development experiments.

They must not be regenerated or overwritten automatically by normal CI.

## 5. Pre-canonical / legacy artifact roots

The repository also retains older `lightning_*` experiment and result roots from the development lineage.

Do not move or delete these directories during package hardening. Their references should be audited first so historical provenance is not broken merely to make the root directory look cleaner.

## 6. Tests

```text
tests/
    test_frozen_v83_parity.py
```

The first hardening gate compares the packaged B2/B3/B4/B5 implementation and router against the historical frozen experiment implementation.

Expected policy:

- no retuning;
- no updated historical expected values;
- no hidden-set optimization;
- strict deterministic parity before further refactoring.

## 7. CI policy

Package-hardening CI should be read-only with respect to historical results.

It may:

- install the package;
- run unit/parity tests;
- run small deterministic smoke tests.

It must not automatically rewrite frozen benchmark outputs.

## 8. Current cleanup order

1. Prove exact V8.3 package parity.
2. Add deterministic unit/integration coverage.
3. Migrate or replace the legacy CLI only after parity is established.
4. Harden reproducibility and dataset provenance.
5. Convert historical workflows to read-only/manual validation.
6. Audit root-level legacy artifacts before any physical directory moves.
7. Only after the engineering preflight is complete should a Phase-2 experimental branch be created.
