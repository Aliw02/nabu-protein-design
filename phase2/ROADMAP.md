# NABU Phase 2 Roadmap

## Execution status

- [x] 2A.0 — Campaign Simulator and Oracle Isolation
- [x] 2A.1 — Baselines
- [x] 2A.2 — Micro-Batch Closed Loop
- [x] 2A.3 — Evidence-Adaptive Controller
- [x] 2A.4 — Multi-Landscape Evaluation
- [x] 2B — Candidate Assembly
- [x] 2C — Assembly + Acquisition (frozen with limitation)
- [ ] 2D — Freeze + Untouched Transfer ← **NEXT**

2A.2 development freeze: use **micro_8** as the primary 2A.3 development cadence and **micro_16** as a sensitivity control. This is a development choice, not a universal optimum claim.

2A.3 development freeze: carry **historical_50_50** and **relative_evidence_v2** into 2A.4 without further RhlA tuning. V1 remains preserved as a failed development attempt.

2A.4 development freeze: multi-landscape acquisition is **PASS WITH LIMITATION**. Gains reproduced on GB1 and TrpB, while PhoQ remains a real failure case. The 2B gate is open only because candidate assembly will be developed separately from acquisition; combined behavior remains a 2C question.

2B development freeze: candidate assembly is **PASS WITH LIMITATION**. The engine is valid and generates novel-to-table combinations, but RhlA quality improves only after evidence matures. 2C must use a generic evidence-state maturity gate; do not hardcode a 20% dataset threshold.

2C development freeze: V2 is **FROZEN WITH LIMITATION**. It passed engineering/leakage checks and increased GB1 Top-1% discoveries from 3 to 10, but did not improve the preregistered best-fitness AUC. No V3 tuning on GB1 is allowed before 2D.

---

## Purpose

Phase 2 changes NABU from static ranking of a supplied candidate library into a controlled closed-loop combinatorial design system.

Phase 1 answered:

> Given measured variants and a candidate pool, can the frozen NABU V8.3 core infer and rank unseen combinatorial variants?

Phase 2 asks:

> Given a known protein, a defined mutation vocabulary, and a limited measurement budget, can NABU iteratively choose or assemble experiments and discover high-fitness variants efficiently?

This file is the execution contract for Phase 2. Future sessions should read this file first before modifying Phase-2 code or running experiments.

---

## Non-negotiable boundary

The frozen Phase-1 V8.3 core remains the reference baseline.

Canonical Phase-1 core:

`freeze/nabu-v8-3-dual-objective-router-validated`

Commit:

`c8afdcd7698231d95a39ef13a3fe22b6e3f507c5`

Phase-1 closure:

`freeze/nabu-phase1-final-cr9114-closed`

Commit:

`690494691c20edc3c9b7a6bba9a5156913173c58`

Phase 2 must not silently retune Phase-1 thresholds, rewrite historical outputs, alter sealed Phase-1 truth, or redefine the meaning of previous Phase-1 results.

Any new mechanism belongs outside the frozen V8.3 core unless it is explicitly versioned as a new architecture.

---

# Execution order

The required order is:

1. 2A.0 — Campaign simulator and oracle isolation
2. 2A.1 — Baselines
3. 2A.2 — Micro-batch closed loop
4. 2A.3 — Evidence-adaptive controller
5. 2A.4 — Multi-landscape evaluation
6. 2B — Candidate assembly
7. 2C — Assembly + acquisition
8. 2D — Freeze + untouched transfer

Do not start candidate assembly before the closed-loop acquisition layer is scientifically understandable.

Do not combine several new mechanisms in one first experiment. If an experiment fails, we must be able to identify which component caused the failure.

---

# 2A — Closed-Loop Acquisition

## 2A.0 — Campaign Simulator and Oracle Isolation

### Goal

Build a true closed-loop virtual wet-lab engine around the frozen V8.3 package.

The candidate pool may be known to the simulator, but the acquisition policy must only see candidate identities and information derivable from already measured data.

The policy must never access unselected truth.

### Required state

A campaign state should track at minimum:

- measured candidate IDs;
- measured labels;
- unmeasured candidate identities;
- mutation sets;
- round number;
- measurements spent;
- current fitted NABU V8.3 model;
- current router mode;
- frozen mutation vocabulary;
- selection hashes / audit metadata.

### Required loop

```text
CampaignState
    ↓
fit frozen V8.3
    ↓
AcquisitionPolicy sees identities + measured evidence only
    ↓
select batch IDs
    ↓
freeze selected IDs
    ↓
VirtualAssayOracle reveals ONLY selected labels
    ↓
update measured set
    ↓
refit
    ↓
next round
```

### Anti-leakage requirements

- Unselected labels are inaccessible to acquisition code.
- Selected candidate IDs are frozen before reveal.
- Truth access is isolated behind a virtual-assay/oracle interface.
- Round manifests record selected IDs and hashes before labels are returned.
- A negative leakage test must deliberately change hidden/unselected truth and prove that pre-reveal selection does not change.

### Bootstrap viability gate

The frozen V8.3 core uses deterministic five-fold cross-fitting and requires all five folds to be populated.

Therefore the initial measured seed cannot simply be an arbitrary tiny random subset.

Before first fit, verify:

- all five deterministic folds are populated;
- labels are finite;
- candidate IDs are stable and unique;
- mutation identities are valid;
- the initial seed is large enough for a valid V8.3 fit.

If the seed fails this gate, extend the bootstrap using identity-only logic. Do not inspect fitness to repair the seed.

### Exit criteria

2A.0 passes only when:

- a complete measure → reveal → update → refit loop runs;
- acquisition code cannot read unselected truth;
- deterministic replay produces the same selections;
- audit manifests reproduce the round sequence;
- the frozen V8.3 package is used without mathematical retuning.

Results belong in:

`phase2/2A_closed_loop_acquisition/2A0_campaign_simulator/results/`

---

## 2A.1 — Baselines

### Goal

Create fair controls before inventing a new controller.

### Required baselines

At minimum:

1. deterministic random acquisition;
2. greedy exploitation;
3. exploration-only;
4. frozen historical 50/50 exploration/exploitation controller;
5. static initial V8.3 ranking without closed-loop adaptation;
6. later, the new Phase-2 adaptive controller.

### Important V8.3 score rule

`V8_3_ADAPTIVE_ROUTER` is a ranking output, not a universally calibrated fitness score.

- `GLOBAL_HIGHER_ORDER`: output uses B5-valued scores.
- `B3_PROTECTED_NO_HIGHER_ORDER`: output uses B3-valued scores.
- `RANK_PRESERVING_B3_TOP20_B5_RERANK`: output is an ordinal rank-preserving score.

Therefore acquisition must not compare raw V8.3 magnitudes across router regimes as though they were on one fitness scale.

For exploitation, use within-pool rank or percentile rank unless a specific calibrated quantity is explicitly justified.

### Fair-comparison rules

All policies must share:

- same initial seed;
- same total measurement budget;
- same candidate pool;
- same oracle;
- same batch schedule;
- same tie-breaking;
- same stopping rules;
- same evaluation checkpoints.

Results belong in:

`phase2/2A_closed_loop_acquisition/2A1_baselines/results/`

---

## 2A.2 — Micro-Batch Closed Loop

### Motivation

The historical 50/50 controller advanced through large budget jumps such as:

`5% → 10% → 20% → 40%`

Within a large jump, many candidates were selected from one model state before refitting.

That creates batch staleness: later selections in the same large batch do not benefit from labels revealed for earlier selections.

This is a plausible explanation for part of the observed non-monotonic active-acquisition behavior, but it is not yet proven.

### Goal

Test whether frequent refitting improves campaign stability and discovery efficiency.

### Required comparison

Compare at least:

- large historical-style batches;
- smaller micro-batches;
- optionally single-candidate rounds where computationally practical.

Do not change the acquisition scoring rule in the first batch-size experiment. First isolate the effect of refit cadence itself.

### Primary question

Does smaller-batch refitting reduce or remove the non-monotonic degradation observed at larger active budgets?

### Exit criteria

We must know whether batch staleness is:

- a major contributor;
- a minor contributor;
- or not supported by the evidence.

Results belong in:

`phase2/2A_closed_loop_acquisition/2A2_micro_batch_loop/results/`

---

## 2A.3 — Evidence-Adaptive Controller

### Goal

Replace a fixed 50/50 controller with a policy whose exploration/exploitation behavior responds to the current state of knowledge.

Do not begin by tuning another fixed scalar mixture such as 0.63/0.37.

### Candidate evidence available from NABU

Potential features include:

- B3 rank;
- B4 rank;
- B5 rank;
- B5 − B3 disagreement;
- main-effect support;
- pair support;
- triplet supported count;
- quartet supported count;
- triplet confidence;
- quartet confidence;
- scoreable / unscoreable status;
- router regime;
- fold-level prediction/rank stability;
- candidate novelty;
- batch diversity.

### Important confidence boundary

Current `triplet_confidence` and `quartet_confidence` are support/coverage confidence measures.

They are not, by themselves, calibrated statistical uncertainty.

Do not label them as epistemic uncertainty.

### Phase-2 uncertainty diagnostic

Use the existing cross-fitting structure to investigate a shadow ensemble diagnostic:

- score a candidate under fold-specific/shadow fitted models;
- measure prediction dispersion;
- measure rank instability;
- keep this outside the frozen V8.3 core.

This can provide a more defensible uncertainty-like signal without modifying Phase-1 mathematics.

### Controller principle

When evidence is weak, support is sparse, or model disagreement is high:

`exploration pressure ↑`

When support is strong and ranking is stable:

`exploitation pressure ↑`

The transition should be evidence-driven, not hardcoded by dataset name.

### Diversity

Acquisition should avoid spending an entire batch on near-duplicate candidates unless the experiment specifically requires it.

Diversity must be treated as a campaign-level selection property, not as hidden biological knowledge.

Results belong in:

`phase2/2A_closed_loop_acquisition/2A3_adaptive_controller/results/`

---

## 2A.4 — Multi-Landscape Evaluation

### Goal

Test the closed-loop controller across multiple landscapes before allowing candidate assembly.

### Required evaluation properties

- fixed preregistered budgets;
- identical campaign protocol across controllers;
- deterministic seeds/tie-breaking;
- no tuning on the held-out outcome being reported;
- per-landscape results plus aggregate summaries;
- no claim that one landscape proves universal active learning.

### Primary campaign metrics

Phase 2 is not primarily a global-correlation benchmark.

Primary metrics should include:

- best true fitness found vs measurements spent;
- measurements to first Top-1% discovery;
- cumulative Top-1% discoveries;
- best-so-far regret;
- area under the discovery curve;
- diversity among high-fitness discoveries;
- measurement efficiency.

Useful diagnostics include:

- Spearman;
- scoreable coverage;
- main/pair/triplet/quartet memory growth;
- router-mode trajectory;
- fold/rank instability trajectory.

### Gate to 2B

Do not begin candidate assembly until we can explain:

- whether closed-loop refitting helps;
- which controller behavior is stable;
- where the controller still fails;
- whether gains reproduce across more than one landscape.

Results belong in:

`phase2/2A_closed_loop_acquisition/2A4_multilandscape_evaluation/results/`

---

# 2B — Candidate Assembly

## Goal

Move from selecting candidates that already exist in a supplied table to assembling novel mutation combinations inside a known mutation vocabulary.

This is still not de novo full-sequence protein generation.

### Conservative initial boundary

- known protein;
- known reference sequence;
- known mutable positions / mutation vocabulary;
- no unknown amino-acid vocabulary;
- no structure model required for the first version;
- no hardcoded dataset-specific biological wordlists;
- adaptive search to avoid combinatorial explosion.

### Required ReferenceProtein layer

Phase-1 token validation can verify that a token has the form `A10C`, but it does not verify that residue 10 in the actual reference protein is `A`.

Phase 2B therefore requires a reference-protein constraint layer that validates substitutions against the true reference sequence.

Example:

If WT position 10 is `A`:

- `A10C` may be valid;
- `G10C` must be rejected for that reference.

### Assembly concept

```text
known mutation vocabulary
        ↓
main / pair / triplet / quartet memory
        ↓
supported partial assemblies
        ↓
adaptive expansion
        ↓
validity + duplication constraints
        ↓
NABU ranking
        ↓
candidate proposals
```

### Replicate handling

Real assays may measure the same variant more than once.

The frozen V8.3 facade intentionally expects unique canonical mutation rows.

Phase 2 therefore needs a separate Measurement Aggregator:

```text
raw assay replicates
        ↓
frozen aggregation rule
        ↓
one canonical candidate label
        ↓
NABU
```

The aggregation rule must be fixed before evaluating a campaign.

Results belong in:

`phase2/2B_candidate_assembly/results/`

---

# 2C — Assembly + Acquisition

## Goal

Combine novel candidate assembly with closed-loop experimental selection.

```text
measured seed
    ↓
fit NABU
    ↓
assemble novel combinations
    ↓
score / filter / diversify
    ↓
choose experimental batch
    ↓
oracle / assay
    ↓
update memory
    ↓
repeat
```

### Main scientific question

Can NABU generate and prioritize useful combinations that were not already supplied as a candidate table, while discovering high-fitness variants under a fixed measurement budget?

### Required controls

Separate the value of:

- assembly alone;
- acquisition alone;
- assembly + acquisition.

Do not attribute improvement to the combined system without these ablations.

Results belong in:

`phase2/2C_assembly_plus_acquisition/results/`

---

# 2D — Freeze + Untouched External Benchmark

## Goal

Freeze the final Phase-2 design/acquisition system and evaluate it on reserved external datasets that have not participated in development.

The authoritative protocol is:

`phase2/FINAL_VALIDATION_PROTOCOL.md`

The machine-readable dataset status is:

`phase2/BENCHMARK_LEDGER.json`

### Development datasets are excluded from final PASS

RhlA, GB1, TrpB, PhoQ, PHOT, eqFP611, CreiLOV and CR9114-H1 remain part of the cumulative evidence record, but they do not count as blind final evidence because they have already influenced development or interpretation.

### Reserved final benchmarks

Primary direct head-to-head:

- FLIP AAV/Random;
- published retrospective Bayesian-optimization comparison against greedy, UCB, Thompson Sampling and random;
- use the same candidate pool, initial folds, measurement checkpoints and oracle;
- primary winner metric: area under Top-100 recovery versus measurements.

Secondary external transfer:

- FLIP2 IRED two-to-many;
- train domain: 0–2 mutations;
- test domain: higher-order variants;
- report Spearman, NDCG, enrichment and regret against published FLIP2 baselines where directly comparable.

### Before any blind target labels are loaded

Freeze:

- 2C architecture;
- maturity gate;
- acquisition policy;
- batch schedule;
- assembly rules;
- reference constraints;
- aggregation policy;
- metrics;
- winner rule;
- success criteria;
- random seeds / folds;
- deterministic tie-breaking;
- source hashes;
- candidate/selection manifests.

Adapters for blind datasets must be debugged on synthetic/development fixtures only.

### After reveal

No scientific retuning is allowed.

Only demonstrably label-independent implementation bug fixes may be made, and each such fix requires:

- preservation of the failed run;
- a new version;
- an audit note;
- complete rerun.

The direct benchmark winner must be reported from the preregistered primary metric. If uncertainty prevents a meaningful distinction, report a tie/unresolved result rather than changing the metric.

Results belong in:

`phase2/2D_untouched_transfer/results/final_blind/`

---

# Metrics contract

The main Phase-2 plot should be conceptually:

```text
Best Fitness Discovered
          ↑
          |
          |
          +--------------------→ Measurements Spent
```

Global Spearman remains useful, but it is a diagnostic rather than the primary objective.

Minimum campaign reporting:

- best-so-far fitness curve;
- first Top-1% discovery budget;
- cumulative Top-1% hits;
- regret curve;
- discovery-curve area;
- diversity;
- scoreability coverage;
- memory/support growth;
- router trajectory.

---

# Engineering / scientific discipline

For every Phase-2 section:

1. preregister the question being tested;
2. isolate one major change whenever possible;
3. save structured JSON diagnostics;
4. save candidate selections before truth reveal;
5. preserve failures, not only successful outputs;
6. do not overwrite prior results;
7. use deterministic tie-breaking;
8. avoid dataset-name hardcoding;
9. do not use hidden truth for acquisition or assembly decisions;
10. keep the frozen Phase-1 core reproducible.

Each section has its own `results/` directory. Never dump all Phase-2 outputs into one shared result folder.

---

# Immediate next step

The current implementation task is **2C — Assembly + Acquisition**.

2C must be completed and frozen using development datasets only.

Do not download or inspect target labels from the reserved AAV/Random or FLIP2 IRED final benchmarks before the 2C freeze commit exists.

After 2C freezes, execute the preregistered 2D blind benchmark exactly as specified in `FINAL_VALIDATION_PROTOCOL.md`.
