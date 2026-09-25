# NABU Repository Map

This document separates the repository by role so the scientific audit trail and the reusable software package are not confused.

## 1. Canonical reusable runtime

```text
src/nabu_protein/
    core.py           Legacy pairwise package API kept for compatibility
    higher_order.py   Frozen Phase-1 B2/B3/B4/B5 hierarchy
    router.py         Frozen V8.3 dual-objective OOF router
    metrics.py        Canonical frozen V8.3 evaluation metrics
    v83.py            Reusable validated facade over the frozen V8.3 core
    v83_cli.py        Explicit convenience CLI for the frozen V8.3 package
    cli.py            Legacy pairwise compatibility CLI
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
    test_v83_contracts.py
    test_metrics.py
    test_historical_stage_a_parity.py
    test_stage_ab_sealing.py
```

The hardening suite covers exact frozen-implementation parity, public API contracts, canonical metrics, sealed Stage-A historical score parity, RhlA budget parity, and Stage-A/Stage-B freeze semantics.

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

## 8. Engineering preflight status

Completed on `engineering/v8-3-package-hardening`:

1. Exact frozen V8.3 implementation parity.
2. Historical parity on PHOT, GB1, TrpB, PhoQ, CreiLOV, eqFP611, RhlA, and CR9114-H1.
3. Deterministic unit and integration coverage.
4. Explicit frozen-V8.3 package facade and convenience CLI.
5. Reproducibility and dataset-provenance manifests.
6. Historical workflows converted to read-only/manual validation.
7. Root-level legacy launcher notes relabelled rather than silently treated as current.

Still intentionally deferred:

- physical relocation/deletion of historical `lightning_*` and frozen result roots;
- any Phase-2 branch or new scientific mechanism;
- merging freeze branches back into `main`.

Physical cleanup remains lower priority than preserving audit provenance.
