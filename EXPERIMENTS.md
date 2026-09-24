# NABU Experiment Lineage — Phase 1 Complete

This file indexes the main experimental progression that produced the frozen Phase-1 architecture.

Historical README, protocol, and RUNBOOK files are intentionally preserved as audit records and should not be rewritten after reveal.

## Early lineage

### V1 — structured memory proof of concept
Initial evidence that discrete mutation memory could improve over simple position means.

### V2 — multi-assay replication
Expanded the approach across multiple assays and confirmed that the signal was not dataset-specific.

### V3 — combinatorial shift
GB1, TrpB, PhoQ, TEV and related blind-design experiments established pairwise residual epistasis as the most reliable backbone.

### V4–V6 — methodology hardening
Introduced stricter candidate freezing, genotype-hash splitting, negative controls, and anti-leakage discipline.

## Modern architecture lineage

### V7 — PHOT public baseline
Higher-order modeling improved global rank but exposed methodological weaknesses in how residual interactions were estimated.

### V8 — true triplet/quartet residual hierarchy
Introduced explicit B2/B3/B4/B5 decomposition with triplet and quartet residual memories.

### V8.1 — strict cross-fitting
Removed same-row self-influence from higher-order residual targets.

On PHOT, higher-order benefit survived strict cross-fitting.

### V8.1 GB1 replication
Higher-order global Spearman did not improve, but elite Top-K discovery did. This showed that global rank and elite discovery are different objectives.

### V8.2 objective tournament
Tested multiple higher-order routing strategies and established that B3 should remain the protected global backbone.

### V8.3 initial adaptive router
Used visible OOF global signal to choose between B3 and higher-order scoring.

Initial PhoQ behavior showed that negative global OOF alone was not enough to decide elite reranking.

### V8.3 Dual-Objective Router
Added an independent visible OOF elite gate.

Final routing modes:

1. `GLOBAL_HIGHER_ORDER`
2. `RANK_PRESERVING_B3_TOP20_B5_RERANK`
3. `B3_PROTECTED_NO_HIGHER_ORDER`

Multi-landscape outcome:

- PHOT -> global HO
- TrpB -> global HO
- GB1 -> elite-only rerank
- PhoQ -> protect B3

This became the canonical Phase-1 core.

## Sealed validation lineage

### CreiLOV sealed extrapolation

Folder: `experiments/v8_3_creilov_sealed`

Purpose: test low-order training against hidden 4/5-mutants.

Outcome:

- strong B3 extrapolation,
- no supported higher-order memory,
- router protected B3,
- Spearman 0.8968.

Scientific label: **sealed pairwise extrapolation validation**.

### eqFP611 sealed higher-order validation

Folder: `experiments/v8_3_eqfp611_sealed`

Purpose: fresh test with supported repeated triplet/quartet contexts.

Stage A:

- 286 triplets
- 715 quartets
- positive visible OOF higher-order signal

Hidden outcome:

- Spearman 0.7601 -> 0.7894
- strict gain on 4/5 primary endpoints
- `SEALED_CANDIDATE_TEST_PASS`

This established fresh higher-order hidden-set benefit.

## Sample-efficiency lineage

### RhlA sealed sample efficiency

Folder: `experiments/v8_3_rhla_sample_efficiency`

Purpose: determine how much measured data is required to approach the 40% reference.

Fixed hidden triple-mutant pool: 540 variants.

Result:

- 5%: Spearman 0.5685
- 10%: 0.6785
- 20%: 0.6900
- 40%: 0.7620

The full frozen practical threshold was first reached at 40%.

### RhlA Active Acquisition development

Folder: `experiments/v8_3_rhla_active_acquisition`

Purpose: test a 50/50 structural-exploration and V8.3-exploitation acquisition controller.

Important boundary: retrospective controller development, not fresh proof.

Result:

- Active 20% reached Passive 40% threshold
- 50% measurement reduction
- 2.0× efficiency factor

The controller was frozen for one final fresh test without retuning.

## Final Phase-1 gate

### CR9114-H1

Folder: `experiments/v8_3_cr9114_final_phase1`

Purpose: combine the frozen V8.3 core and frozen Active Acquisition controller in one fresh high-order 4/5-mutation gate.

Before reveal:

- hidden 4/5-mutants: 1,673
- coverage: 100%
- Passive 40 triplets/quartets: 546 / 1,698
- Active 20 triplets/quartets: 545 / 1,024

Core outcome:

- B3 Spearman 0.9271
- V8.3 Spearman 0.9387
- 4/5 strict gains
- zero regression
- **core PASS**

Active outcome:

- Active 10% reached Passive 40% threshold -> 4× efficiency
- Active 20% failed the preregistered final target
- active behavior was non-monotonic

Formal status:

`PHASE1_FINAL_GATE_CLOSED_WITH_LIMITATION`

No post-reveal retuning is allowed.

## Frozen Phase-1 references

Core architecture:

`freeze/nabu-v8-3-dual-objective-router-validated`  
`c8afdcd7698231d95a39ef13a3fe22b6e3f507c5`

Final closure:

`freeze/nabu-phase1-final-cr9114-closed`  
`690494691c20edc3c9b7a6bba9a5156913173c58`

## Next

Phase 2 begins with closed-loop combinatorial design rather than another static ranking benchmark.

See [PHASE2_PLAN.md](PHASE2_PLAN.md).
