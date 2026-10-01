# NABU Scientific Constraints and Canonical Reset State

**Canonical branch:** `main`  
**Reset date:** 2026-10-01  
**Status:** Phase 3 CLOSED · NucB UNTOUCHED · replacement architecture NOT FROZEN

This file is the authoritative scientific constraint document for all work after the Phase-2D transfer failure. If another experiment note, branch README, or prior chat conflicts with this file, this file governs until it is explicitly revised on `main`.

## 1. Canonical rollback checkpoint

Development is rolled back conceptually to the last defensible baseline plus the root-cause falsification evidence.

### Frozen baseline

TrpB blind B2 remains the reference transfer baseline:

- full higher-order Spearman: **0.4706975858**
- triples Spearman: **0.6751217397**
- quads Spearman: **0.1778461371**

Interpretation: the system retained useful higher-order ranking signal, but generalization degraded sharply from triples to quads.

### Root-cause status

**ROOT_CAUSE_NOT_YET_IDENTIFIED**

Do not promote a diagnostic observation into a mechanistic root cause unless a preregistered discriminative test separates it from plausible competing explanations on more than one revealed development landscape.

## 2. What is scientifically supported

The following observations are retained as canonical evidence.

### Implementation/data-semantics failure is not the main explanation

Stage-1 invariants passed on revealed TrpB and IRED:
- exact low-order reconstruction,
- mutation parsing/reference handling,
- support bookkeeping,
- metric self-tests,
- split/source identity checks.

### Basic low-order directional representation failure is not the main explanation

Observed alignment between exact low-order quantities and learned low-order terms was material:

- TrpB singleton vs B2 main Spearman: **0.7998737410**
- TrpB pair e2 vs B3 pair contribution: **0.7093196656**
- IRED singleton vs B2 main: **0.6707548**
- IRED pair e2 vs B3 pair contribution: **0.7251674**

Scale shrinkage exists, but the directional structure is not absent.

### Context-free pair portability is not supported

On both revealed TrpB and GB1, an exact pair effect measured in a double could not be treated as a fixed additive contribution inside higher-order variants.

### Conditional interaction context dependence is material

Observed-scale diagnostics on exact-O2-supported triples showed:

- TrpB pair sign-flip rate: **0.3630964433**
- TrpB e3 RMS / target SD: **5.9524430494**
- GB1 pair sign-flip rate: **0.2695160109**
- GB1 e3 RMS / target SD: **1.6838263644**

### Observation-scale nonlinearity contributes but is not sufficient

A low-order-fit robust-asinh transform materially improved some triple reconstruction metrics, but did not remove pair non-portability or context dependence on both landscapes.

### A single scalar cancellation law is not universal

A source-triple affine cancellation law improved TrpB quads but did not transfer equivalently on GB1.

### Material fourth-order residual remains

After scale correction and measured O3 reconstruction:

- TrpB e4 RMS / target SD: **2.8157971835**
- GB1 e4 RMS / target SD: **1.6952732515**

### Identifiability boundary

WT + singles + doubles do not uniquely identify unrestricted third-order or higher Möbius coefficients without an additional structural prior.

This is a mathematical information constraint. It does **not** imply that useful zero-shot higher-order prediction is impossible; it means such prediction requires justified prior structure or additional measurements.

## 3. What is NOT established

The following claims must not be stated as facts:

- that the biological root cause has been identified;
- that low-order labels alone are sufficient to reconstruct arbitrary higher-order interactions;
- that scale nonlinearity alone explains the transfer failure;
- that interaction-shape statistics are a validated replacement architecture;
- that a random-forest residual transfer result proves target low-order sufficiency;
- that a source-protein higher-order prior solves the original target-only inference problem.

## 4. Retired side investigations

The following results remain audit evidence but are **not** the active architecture path:

- cross-landscape random-forest higher-order prior transfer;
- cross-order prior transfer;
- final interaction-shape ablation;
- NABU V10.0 B2 + external interaction-shape rank correction.

V10.0 was explicitly rejected by its frozen development gate.

Canonical V10.0 rejection:
- run: `36418960080`
- artifact: `10968927659`
- result commit: `be9aa004423a079906e7d59a33d0d87ff9102d25`

Do not create V10.1 as a repair of the same rank-correction architecture.

## 5. Inference-time information contract

For a new target landscape, every experiment must state exactly what information is available at inference time.

Unless a future preregistered protocol explicitly changes the scientific question, the target-side supervised evidence is limited to:

- WT,
- measured singles,
- measured doubles,
- mutation vocabulary / sequence identity information that is available without using hidden target fitness labels.

Target triple/quad/higher-order fitness labels are forbidden for:
- fitting,
- feature selection,
- transform selection,
- architecture selection,
- threshold selection,
- early stopping,
- calibration,
- blend-weight selection.

They may be used only once for revealed-development evaluation when the protocol explicitly designates that dataset as revealed development.

## 6. NucB protection

NucB remains the untouched external benchmark.

Before NucB can be consumed, all of the following must already be frozen:

1. architecture;
2. feature set;
3. training procedure;
4. external/pretrained prior, if any;
5. hyperparameters;
6. support/fallback policy;
7. evaluation metrics;
8. PASS/FAIL gate;
9. source hashes;
10. deterministic replay check.

No NucB result may cause any of those items to change.

If a NucB run fails, the failure is recorded. NucB is not converted into a development dataset.

## 7. Forbidden development behavior

Do not:

- patch a failed architecture with post-reveal weights or thresholds;
- create sequential V10.1/V10.2-style repairs from the same revealed failure;
- add arbitrary blend coefficients;
- add order-specific corrections solely because order-3 or order-4 failed;
- add a local reranker and call it a root-cause fix;
- treat exact e3/e4 measured on revealed targets as information available at inference time;
- mix target hidden labels into feature normalization or transform fitting;
- change a preregistered gate after seeing the result;
- delete failed results;
- describe GitHub workflow success as scientific PASS;
- call an observation a root cause unless competing hypotheses were discriminated prospectively;
- open Phase 3 because one development metric improved.

## 8. Required experiment admission checklist

No new architecture experiment should run unless its protocol, committed before execution, contains all of the following:

1. **Scientific question** — one sentence.
2. **Competing hypotheses** — including what result would reject each.
3. **Inference-time information** — exact allowed inputs.
4. **Forbidden information** — exact leakage boundary.
5. **Baseline** — frozen comparator.
6. **Intervention** — one architectural change only.
7. **Controls** — at least one null/shortcut control when applicable.
8. **Dataset roles** — development vs untouched external.
9. **Metrics** — declared before reveal.
10. **PASS/FAIL gate** — declared before reveal.
11. **Kill rule** — what is permanently abandoned on failure.
12. **Determinism** — seed/replay requirement.
13. **Artifacts** — JSON/CSV outputs, manifests, source hashes.
14. **No-patch rule** — explicit statement that failure will not be repaired on the same revealed target.

If any item is missing, the experiment is not admitted.

## 9. Branch discipline

Use `main` for canonical frozen state and scientific constraints.

At most one active architecture-development branch should exist at a time.

Before deleting an experimental branch:
- preserve the canonical result summary and immutable commit/run/artifact identifiers on `main`;
- classify it as PASS, REJECT, or NON-ADJUDICATING;
- preserve negative results.

Branches designated for cleanup after this reset:

- `exp/nabu-root-cause-audit` — retired after its canonical evidence is captured here;
- `exp/nabu-v10-interaction-shape-dev` — retired; V10.0 rejected.

Historical Phase-1/Phase-2 branches are audit history and are not part of this cleanup unless separately reviewed.

## 10. Current scientific question

The project is not currently searching for another scalar correction.

The unresolved question is:

> What justified structural prior, available without target higher-order fitness labels, can explain or predict the context-dependent change in interaction effects strongly enough to improve higher-order transfer beyond the frozen B2 baseline?

A sequence-aware or pretrained context prior is one candidate class, not a predetermined answer.

No implementation of that class should begin until a protocol satisfies the experiment admission checklist above.

## 11. Current state

- Phase 1: CLOSED and preserved.
- Phase 2A/2B/2C: preserved as completed historical work.
- Phase 2D external transfer: scientific replacement gate NOT PASSED.
- B2: frozen transfer baseline.
- Root cause: NOT YET IDENTIFIED.
- V10.0: REJECTED.
- NucB: UNTOUCHED.
- Phase 3: CLOSED.

The next action is scientific problem formulation under this document, not another architecture patch.
