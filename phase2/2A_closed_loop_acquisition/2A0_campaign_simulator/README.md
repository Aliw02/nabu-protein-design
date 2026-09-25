# 2A.0 — Campaign Simulator and Oracle Isolation

## Status

Implementation scaffold for the first Phase-2 gate.

This stage does **not** introduce a new active-learning controller and does not make a new scientific performance claim.

Its purpose is to prove that NABU can run in a deterministic, replayable, leak-resistant closed loop around the frozen V8.3 core.

## Files

- `campaign.py` — generic campaign state, identity-only bootstrap, policy view, frozen selection, oracle reveal, refit, and audit manifests.
- `run_rhla_smoke.py` — engineering smoke run using the RhlA training pool as a virtual assay pool.
- `tests/test_campaign_simulator.py` — bootstrap, truth-isolation, pre-reveal freeze, hidden-truth perturbation, and deterministic replay tests.
- `results/` — stage-specific output only.

## Information boundary

The acquisition policy receives a `PolicyView` containing:

- measured candidate IDs;
- current round and budget;
- current router mode;
- unmeasured candidate identities;
- model-derived candidate scores/support diagnostics.

It does **not** receive unmeasured labels.

Only `VirtualAssayOracle.reveal(selected_ids)` returns assay values, after the selected ID list has been frozen and hashed.

## Bootstrap

The initial seed is identity-only.

Candidates are chosen deterministically from candidate IDs while forcing representation in all five frozen cross-fit folds.

No fitness values are used to repair fold coverage.

## Audit files

A run writes:

- `bootstrap/selection_before_reveal.json`
- `bootstrap/manifest.json`
- `rounds/round_XXX_selection_before_reveal.json`
- `rounds/round_XXX_manifest.json`
- `CAMPAIGN_MANIFEST.json`

Selection hashes are created before label reveal.

## RhlA smoke boundary

The RhlA smoke uses only `rhla_sample_efficiency_sealed_input/TRAINING_POOL.csv` as a virtual campaign pool.

The script separates identity columns from labels in memory and gives labels only to the oracle.

It does not use the frozen RhlA hidden-triple truth and should not be interpreted as a new scientific benchmark.

## Pass criteria

2A.0 passes when:

1. all five folds are populated by identity-only bootstrap;
2. policy views contain no hidden truth columns;
3. selections are frozen before reveal;
4. changing unrevealed truth cannot change pre-reveal selection;
5. deterministic replay reproduces the same selection sequence;
6. the frozen V8.3 core refits successfully after each reveal.
