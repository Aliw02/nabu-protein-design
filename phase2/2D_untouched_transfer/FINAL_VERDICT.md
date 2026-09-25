# NABU Phase 2D — Untouched External Benchmark Verdict

## Status

**PHASE 2 FINAL SCIENTIFIC GATE NOT PASSED**

Phase 2D execution is complete under the frozen Phase-2C V3 architecture, with
one label-blind external-scope failure and one completed negative blind transfer
result.

This does not invalidate the Phase-2C development PASS. It limits the external
generalization claim that Phase 2 can support.

## 2D-A — FLIP AAV/Random

Status:

**ADAPTER_SCOPE_FAIL — LABEL BLIND**

The pinned published Bayesian-optimization pool contains 59,459 identities.

Frozen NABU V8.3 can represent:

- 27,018 / 59,459 = 45.44%.

It cannot represent:

- 32,441 / 59,459 = 54.56%.

The unsupported identities include length-changing variants. Frozen Phase-2
mutation semantics do not define insertion identities.

The preregistered protocol forbids silently reducing the published candidate
pool. Therefore:

- no reduced-subset head-to-head was run;
- AAV target labels were not revealed;
- the required primary AAV Top-100-recovery benchmark was not completed.

This is an external benchmark scope limitation, not an AAV fitness-performance
loss.

## 2D-B — FLIP2 IRED two-to-many

Status:

**BLIND RESULT COMPLETE — NEGATIVE HIGHER-ORDER TRANSFER**

Authoritative execution:

- workflow run: `36165359008`;
- execution head: `5c4b55e60bb0ccca769194e78095051e7023b9d2`;
- artifact ID: `10877951045`;
- artifact digest:
  `sha256:32e9b34d6bc7453a54ac1699c5d22a8516185d925a822152f2538770bb22a461`;
- source SHA256:
  `aa45a2f85fb1af87b6b0e86397b3f8f292a061dc574be2536673b5bd490b0e74`.

Frozen split:

- fit/train: 3,746;
- validation held out: 662;
- test: 4,178.

### Label-blind scoreability result

Before target reveal:

- main identities in fit evidence: 1,412;
- pair identities in fit evidence: 2,708;
- full test scoreable by frozen V8.3: **131 / 4,178 = 3.135%**;
- unscoreable: **4,047 / 4,178**;
- test rows with all main mutations supported: 3,648;
- test rows with at least one supported internal pair: only 133.

The dominant bottleneck is therefore pair-combination transfer, not merely
unseen single-mutation identities.

### Authoritative full-test policy

The target-blind policy frozen before reveal was `IRED_ABSTENTION_V2`:

- frozen-scoreable rows use unchanged `V8_3_ADAPTIVE_ROUTER`;
- frozen-unscoreable rows remain abstentions and tie below all scoreable rows;
- no B2 rescue/fallback;
- no test row is deleted;
- primary metrics use the complete 4,178-row test.

### Blind result

Full-test primary metrics:

- Spearman: **-0.0217705732**
- NDCG: **0.9426983698**

Secondary metrics:

- predicted Top-1% budget: 42 variants;
- true Top-1% hits: **0 / 42**;
- Top-1% recall: **0.0**;
- Top-1% enrichment: **0.0**;
- normalized regret at Top-1% budget: **0.2315223517**.

Diagnostic only, on the 131 frozen-scoreable rows:

- Spearman: 0.1193882453;
- NDCG: 0.8771807667.

These diagnostic values do not replace the preregistered full-test result.

### Published FLIP2 context

FLIP2 Table A5 reports, for the same IRED two-to-many split:

- Ridge one-hot: Spearman 0.193, NDCG 0.960;
- Ridge one-hot + likelihoods: Spearman 0.211, NDCG 0.964;
- Dayhoff likelihood: Spearman 0.142, NDCG 0.955;
- ESM2-650M likelihood: Spearman 0.141, NDCG 0.953;
- CARP-640M likelihood: Spearman 0.144, NDCG 0.956.

The untouched frozen NABU result is below these directly comparable published
ranking baselines.

## What 2D establishes

Engineering integrity:

**PASS**

- Phase-2C was frozen before blind contact;
- source hashes were pinned;
- label-blind capability checks preceded target reveal;
- failed contract/preflight runs were preserved;
- no post-reveal scientific retuning was performed.

Active-optimization external validation:

**UNRESOLVED / NOT EXECUTABLE UNDER FROZEN AAV SEMANTICS**

Higher-order external transfer:

**NEGATIVE RESULT**

Assembly evidence:

**2C DEVELOPMENT PASS remains valid, but is not blind final evidence.**

## Final Phase-2 verdict

The preregistered final Phase-2 PASS required both the primary AAV external
head-to-head and the IRED higher-order transfer evaluation.

That condition is not satisfied.

> **NABU Phase 2: FINAL SCIENTIFIC GATE NOT PASSED.**

The correct next research action is not to retune Phase 2 against these blind
results. Any architecture that adds insertion semantics, changes scoreability,
or introduces a new evidence-backoff mechanism must be versioned as new work
outside the frozen Phase-2 claim.
