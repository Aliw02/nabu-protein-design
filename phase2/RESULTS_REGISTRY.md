# Phase 2 Results Registry

This file is the human-readable index for Phase-2 evidence.

The machine-readable benchmark classification is in `BENCHMARK_LEDGER.json`.

## Evidence hierarchy

### Development evidence

- 2A.0: leak-safe campaign simulator.
- 2A.1: acquisition baselines on RhlA.
- 2A.2: micro-batch cadence ablation on RhlA.
- 2A.3: adaptive-controller V1/V2 development on RhlA.
- 2A.4: multi-landscape acquisition evaluation on GB1, TrpB and PhoQ.
- 2B: reference-aware candidate assembly and evidence-maturity analysis on RhlA.
- 2C V3: dynamic maturity-gated assembly + acquisition benchmark on GB1, TrpB and PhoQ; development verdict PASS.

These results remain part of the total scientific record but cannot be used as final blind evidence.

## Mandatory artifact policy

Every scientific run must retain, at minimum:

- source dataset identity and SHA256;
- repository commit SHA;
- frozen configuration;
- seed / fold IDs;
- candidate-selection or proposal hashes before truth reveal;
- full primary metric table;
- summary manifest;
- failed-run status when applicable.

Large per-round logs may be stored as CI artifacts rather than Git, but the artifact ID/digest and summary manifest must be preserved.

## No selective deletion

Negative and neutral results are retained.

Examples already retained:

- 2A.3 V1 underperformance;
- 2A.4 PhoQ failure case;
- 2B weak 5% and 10% assembly quality.
- 2C multilandscape V1 permanent-latch engineering failure.
- earlier GB1-only 2C V1/V2 attempts, retained but superseded for current 2C status by multilandscape V3.

## Final blind evidence

Phase 2D has been executed under the frozen Phase-2C V3 protocol.

### AAV/Random

- status: **ADAPTER_SCOPE_FAIL — LABEL BLIND**;
- identity-preflight workflow: `36164741113`;
- identity-preflight artifact: `10876827869`;
- artifact digest: `sha256:c6ed4dcd89258edf37886d21a2cdb95964a8457a88c3a0eefa1c668fb1279132`;
- published BO pool: 59,459 identities;
- representable by frozen Phase-2 semantics: 27,018 (45.44%);
- unrepresentable: 32,441;
- target labels revealed: **no**;
- primary published AAV head-to-head executed: **no**.

The preregistered no-subsetting rule was preserved.

### FLIP2 IRED two-to-many

- status: **BLIND RESULT COMPLETE — NEGATIVE HIGHER-ORDER TRANSFER**;
- workflow run: `36165359008`;
- execution head: `5c4b55e60bb0ccca769194e78095051e7023b9d2`;
- artifact ID: `10877951045`;
- artifact digest: `sha256:32e9b34d6bc7453a54ac1699c5d22a8516185d925a822152f2538770bb22a461`;
- source SHA256: `aa45a2f85fb1af87b6b0e86397b3f8f292a061dc574be2536673b5bd490b0e74`;
- full test: 4,178 variants;
- frozen-V8.3 scoreable: 131 (3.135%);
- full-test Spearman: **-0.0217705732**;
- full-test NDCG: **0.9426983698**;
- Top-1% discoveries: **0 / 42**;
- normalized regret at Top-1% budget: **0.2315223517**.

Compact authoritative outputs:

- `phase2/2D_untouched_transfer/results/final_blind/IRED_BLIND_V2_SUMMARY.json`
- `phase2/2D_untouched_transfer/results/final_blind/PHASE2D_FINAL_VERDICT.json`
- `phase2/2D_untouched_transfer/FINAL_VERDICT.md`

Full predictions and blind-run diagnostics remain preserved in CI artifact
`10877951045`.

### Final Phase-2 scientific status

Engineering integrity: **PASS**.

Final scientific gate: **NOT PASSED**.

The 2C V3 development PASS remains valid, but the preregistered final Phase-2
external-validation requirements were not satisfied. No Phase-2 scientific
retuning is permitted after these blind outcomes.

## Phase 2C V3 authoritative record

- workflow run: `36160574928`
- validated head: `316910cfc96a4bb98027b01c6f7cbba91646edc2`
- artifact ID: `10876390166`
- artifact digest: `sha256:d12852968cd865775e34abd3474f04e4226e98440b70b6093e8f48763b5c152d`
- deterministic replay: PASS on GB1 / TrpB / PhoQ
- hidden-truth first-decision perturbation: PASS on GB1 / TrpB / PhoQ
- overall development verdict: **PASS**

Canonical compact results:

`phase2/2C_assembly_plus_acquisition/results/multilandscape_v3/`

Full per-round outputs remain preserved in the CI artifact.
