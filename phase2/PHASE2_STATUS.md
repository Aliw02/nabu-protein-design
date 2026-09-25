# NABU Phase 2 Status

## Current state

Completed engineering/scientific development stages:

- **2A.0 — Campaign Simulator and Oracle Isolation:** PASS
- **2A.1 — Baselines:** PASS
- **2A.2 — Micro-Batch Refit Ablation:** PASS
- **2A.3 — Evidence-Adaptive Controller:** DEVELOPMENT COMPLETE
- **2A.4 — Multi-Landscape Evaluation:** PASS WITH LIMITATION
- **2B — Candidate Assembly:** PASS WITH LIMITATION
- **2C — Assembly + Acquisition:** PASS / V3 FROZEN

Next stage:

- **2D — Untouched External Benchmark**

## Frozen Phase-1 dependency

Frozen V8.3 scientific core:

`c8afdcd7698231d95a39ef13a3fe22b6e3f507c5`

Phase-1 closure:

`690494691c20edc3c9b7a6bba9a5156913173c58`

Phase 2 must continue to treat the Phase-1 core as frozen.

## 2A.0 result

The leak-safe campaign engine passed:

- identity-only bootstrap viability;
- five-fold population gate;
- policy/oracle separation;
- pre-reveal selection hashing;
- unrevealed-truth perturbation test;
- deterministic replay;
- reveal → update → refit loop.

## 2A.1 result

RhlA development-landscape controls established:

- deterministic random;
- greedy current V8.3;
- exploration-only;
- historical 50/50;
- static initial V8.3.

Important interpretation:

- discovery quality is not equivalent to global Spearman;
- static ranking underperformed adaptive closed-loop acquisition;
- exploration and exploitation showed complementary behavior;
- the historical 50/50 mixture was a strong development baseline on RhlA.

## 2A.2 result

The acquisition policy was held fixed as `historical_50_50`.

Only reveal/refit cadence changed.

Observed RhlA results:

- historical large jumps: strong degradation;
- micro_32: recovered the global best but late;
- micro_16: strong discovery efficiency;
- micro_8: strongest observed RhlA development result across the preregistered primary metrics;
- micro_4: degraded relative to micro_8 and micro_16.

Conclusion:

> Batch staleness is supported on RhlA. More frequent refitting is useful, but the benefit is not monotonic as batch size approaches very small values.

## Frozen development decision for 2A.3

For **2A.3 controller development only**:

- primary development cadence: **micro_8**;
- sensitivity control: **micro_16**;
- historical 50/50 remains the fixed reference baseline;
- no new fixed scalar mixture may be tuned as the main controller;
- the first 2A.3 experiment must change the controller logic only;
- batch cadence must remain fixed within each comparison;
- no dataset-name-specific rules;
- no hidden truth access;
- no Phase-1 retuning.

This is a development decision, not a universal claim that batch size 8 is optimal.

## Gate before 2B

2A gate decision:

- **OPEN TO 2B WITH LIMITATION**
- active acquisition gains reproduced on GB1 and TrpB;
- PhoQ remains an explicit acquisition failure case;
- 2B must isolate candidate assembly from acquisition;
- assembly + acquisition remains deferred to 2C.


## 2A.3 result

Two adaptive-controller versions were evaluated on RhlA.

V1:
- engineering PASS;
- scientific development underperformance versus historical 50/50;
- failure preserved and diagnosed.

V2:
- self-referenced evidence-gap calibration;
- starts at exploration pressure 0.5;
- pressure falls as evidence gaps shrink;
- deterministic and dataset-name independent;
- improved over V1;
- did not uniformly beat historical 50/50 on RhlA.

Freeze decision:
- carry `historical_50_50` as the fixed reference;
- carry `relative_evidence_v2` as the adaptive candidate;
- do not create V3 by further tuning on RhlA before 2A.4;
- primary cadence remains micro_8;
- micro_16 remains a sensitivity control.


## 2A.4 result

Pinned multi-landscape evaluation completed on GB1, TrpB and PhoQ with no per-landscape retuning.

- GB1: historical 50/50 beat random at micro_8; V2 did not improve the primary micro_8 result.
- TrpB: historical 50/50 and V2 both beat random and tied each other.
- PhoQ: random acquisition outperformed both active controllers at micro_8; V2 improved over historical at micro_16 but did not remove the landscape failure.

Conclusion:

> Active acquisition gains reproduce on more than one landscape, but controller behavior remains landscape-sensitive.

2B Candidate Assembly is allowed to begin as a separate mechanism. The PhoQ acquisition limitation must remain explicit until 2C/2D.


## 2B result

Candidate Assembly V1 is reference-aware, deterministic and leak-resistant.

RhlA assembly POC:
- 64 proposals;
- 64 scoreable;
- 0 measured overlap;
- 0 duplicates;
- hidden-truth perturbation invariant PASS;
- 38/64 proposals were novel to the historical candidate table.

Evidence-maturity result:
- 5% evidence: 0 matched Top-1% hits;
- 10% evidence: 0 matched Top-1% hits;
- 20% evidence: 2 matched Top-1% hits and best known fitness 0.308292.

The 20% improvement occurred with B3-protected routing and no triplet/quartet memories, so denser pairwise evidence was sufficient to improve assembly ranking.

2C must use a generic evidence-state maturity gate before enabling assembly. A fixed 20% threshold is not frozen as universal.


## 2C result

The authoritative Phase-2C result is **V3 Dynamic Maturity Gate**.

Development suite:

- GB1;
- TrpB;
- PhoQ.

The exact same 2048-reservoir / 64-bootstrap / 256-budget / batch-8 protocol was applied without per-landscape tuning.

V1 multilandscape failed engineering because a permanently latched maturity gate could not guarantee a full assembly batch on later rounds. The failure is preserved.

V3 changed only gate persistence:

- assembly is used on a round only when the current mature frontier can supply the complete batch;
- otherwise that round falls back to the frozen acquisition policy.

V3 integrity:

- 17 contract tests PASS;
- full deterministic replay PASS on all 3 landscapes;
- first-decision hidden-truth perturbation invariant PASS on all 3 landscapes;
- no reserved 2D blind target labels loaded.

Primary combined system `acquisition_then_gated_assembly` beat `acquisition_only` under the preregistered rule on all three landscapes:

- GB1: ΔAUC +0.071326, Δregret -1.999854, ΔTop-1% hits +6;
- TrpB: ΔAUC +0.038913, Δregret -0.065887, ΔTop-1% hits +10;
- PhoQ: ΔAUC +0.026854, Δregret -4.265920, ΔTop-1% hits +4.

Overall 2C development verdict:

> **PASS**

V3 is now frozen. No additional tuning on GB1, TrpB or PhoQ is allowed before 2D.
