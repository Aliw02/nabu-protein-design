# Phase 2C V3 — Dynamic Maturity Gate Multilandscape Preregistration

## Lineage

- V1 multilandscape run failed because the assembly gate latched permanently after its first ready round.
- An older GB1-only V2 experiment exists in repository history and changes maturity probe count. It is preserved and is not overwritten.
- V3 addresses only the multilandscape V1 persistence failure.

## Single changed rule versus the failed multilandscape V1

Failed V1:

- once the maturity gate opened, assembly remained permanently active.

V3:

- assembly readiness is evaluated independently on every round;
- assembly is active only when the current mature proposal count is at least the requested experimental batch size;
- if the current mature proposal count is insufficient, that round falls back to the condition's original reservoir acquisition policy.

Therefore:

### random_then_gated_assembly

- mature batch available -> assembly;
- otherwise -> deterministic random reservoir acquisition.

### acquisition_then_gated_assembly

- mature batch available -> assembly;
- otherwise -> historical 50/50 reservoir acquisition.

## Frozen unchanged settings

- development landscapes: GB1, TrpB, PhoQ;
- biological reference genotypes: VDGV, VFVS, AVST;
- supplied reservoir size: 2048;
- bootstrap: 64;
- total measurement budget: 256;
- batch size: 8;
- assembly target order: 3;
- beam width: 256;
- retained proposal frontier: 1024;
- main support >= 2;
- every internal pair support >= 2;
- outside-reservoir requirement;
- assay-universe identity coverage requirement;
- all four conditions;
- all primary metrics;
- per-landscape scientific win criterion;
- overall PASS / PASS_WITH_LIMITATION / NO_SCIENTIFIC_WIN rules;
- deterministic replay requirement;
- state-matched first-ready assembly versus acquisition diagnostic.

No support threshold, search budget, metric, dataset or scientific success rule is changed.

## Dynamic-gate diagnostics

For each assembly-capable condition, record:

- first ready measurement count;
- assembly-ready rounds;
- assembly-not-ready rounds;
- fallback-to-acquisition rounds after first readiness;
- ready/not-ready transition count;
- mature proposal count per round.

## Interpretation

Dynamic gating is a safety condition, not a tuned acquisition hyperparameter:

> use assembly only when the frozen maturity rule can currently supply the full preregistered batch.

V3 outcomes must be preserved even if they fail scientifically.
