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
- 2B: reference-aware candidate assembly and evidence-maturity analysis on RhlA.\n- 2C V3: dynamic maturity-gated assembly + acquisition benchmark on GB1, TrpB and PhoQ; development verdict PASS.

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
- 2B weak 5% and 10% assembly quality.\n- 2C multilandscape V1 permanent-latch engineering failure.\n- earlier GB1-only 2C V1/V2 attempts, retained but superseded for current 2C status by multilandscape V3.

## Final blind evidence

Reserved and not yet loaded:

1. FLIP AAV/Random — direct published active-optimization head-to-head.
2. FLIP2 IRED two-to-many — higher-order transfer challenge.

The final blind outputs will be stored under:

`phase2/2D_untouched_transfer/results/final_blind/`

and will not replace any development results.


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
