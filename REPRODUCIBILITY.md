# NABU Reproducibility Policy

## Scope

Phase-1 historical result files are immutable audit artifacts. Package hardening must not rewrite them.

The repository did not lock one exact dependency environment for every historical Phase-1 experiment. Therefore this document does **not** claim to reconstruct an unknown historical environment.

## Current package-hardening baseline

The engineering hardening gate uses Python 3.11 and the exact direct/test dependencies listed in `requirements-hardening.txt`.

The first successful package-parity run used:

```text
Python 3.11.16
numpy 2.4.6
pandas 3.0.6
scipy 1.17.1
pytest 9.1.1
```

Under that environment, the packaged B2/B3/B4/B5 hierarchy and V8.3 routing behavior matched the frozen experiment implementation.

## Scientific source of truth

Canonical frozen Phase-1 core:

```text
freeze/nabu-v8-3-dual-objective-router-validated
c8afdcd7698231d95a39ef13a3fe22b6e3f507c5
```

Historical experiment code under `experiments/` remains an audit record and parity oracle. It is not silently rewritten during package cleanup.

## Dataset provenance

See `DATASET_PROVENANCE.json`.

Sealed-input manifests remain authoritative for source-file hashes, split definitions, and Stage-A/Stage-B integrity hashes.

## CI policy

Normal engineering CI may:

- install the package;
- run deterministic unit, integration, and parity tests;
- create ephemeral CI artifacts.

Normal engineering CI must not:

- rewrite historical benchmark result directories;
- mutate preregistration files;
- mutate sealed-input manifests;
- retune frozen scientific thresholds.
