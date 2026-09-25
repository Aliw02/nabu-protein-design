# 2D — Untouched External Benchmark

## Final status

**COMPLETE — PHASE-2 FINAL SCIENTIFIC GATE NOT PASSED**

Phase 2D was executed only after the Phase-2C V3 freeze.

Authoritative Phase-2D verdict:

`FINAL_VERDICT.md`

Machine-readable verdict:

`results/final_blind/PHASE2D_FINAL_VERDICT.json`

## AAV/Random

Status:

`ADAPTER_SCOPE_FAIL — LABEL BLIND`

The published BO pool contains 59,459 identities. Frozen Phase-2 mutation
semantics represent 27,018 (45.44%) and cannot represent 32,441 length-changing
variants without adding new insertion semantics.

The preregistered no-subsetting rule was preserved. AAV target labels were not
revealed and no reduced-pool head-to-head was reported.

## FLIP2 IRED two-to-many

Status:

`BLIND RESULT COMPLETE — NEGATIVE HIGHER-ORDER TRANSFER`

Authoritative run:

- workflow: `36165359008`
- execution head: `5c4b55e60bb0ccca769194e78095051e7023b9d2`
- artifact: `10877951045`
- artifact digest: `sha256:32e9b34d6bc7453a54ac1699c5d22a8516185d925a822152f2538770bb22a461`

Full-test result:

- scoreable coverage: 131 / 4,178 (3.135%)
- Spearman: -0.0217705732
- NDCG: 0.9426983698
- Top-1% hits: 0 / 42
- normalized regret: 0.2315223517

## Scientific boundary

The Phase-2C V3 development PASS remains valid.

The final Phase-2 external scientific gate was not passed.

Do not retune the frozen Phase-2 system against these blind outcomes. Any
insertion support, changed scoreability rule, or new evidence-backoff mechanism
must be versioned as new research.
