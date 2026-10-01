# NABU Scientific Evidence Ledger

This file is the canonical index of post-reset scientific evidence. Detailed experiment records live under `evidence/`. Negative results are permanent evidence, not disposable failures.

## Canonical baseline

### TrpB blind B2

- full higher-order Spearman: 0.4706975858
- triples Spearman: 0.6751217397
- quads Spearman: 0.1778461371

Status:
FROZEN BASELINE

## Root-cause falsification checkpoint

Status:
ROOT_CAUSE_NOT_YET_IDENTIFIED

Retained conclusions:
- implementation/data semantics are not the main explanation;
- basic low-order directional representation is materially present;
- context-free pair portability is not supported on TrpB and GB1;
- conditional interaction context dependence is material;
- robust-asinh scale correction contributes but is insufficient;
- one scalar cancellation law is not universal;
- material fourth-order residual remains;
- WT + singles + doubles do not uniquely identify unrestricted third-order or higher terms without an additional prior.

## V10.0 interaction-shape rank correction

Run:
`36418960080`

Artifact:
`10968927659`

Decision:
`V10_DEV_REJECT`

Meaning:
Transferable interaction-shape structure was measurable in revealed datasets, but the concrete B2 + external interaction-shape global rank-correction architecture did not transfer robustly.

Status:
RETIRED ARCHITECTURE

## Sequence-Context Falsification 01

Detailed canonical record:
[evidence/SEQUENCE_CONTEXT_FALSIFICATION_01.md](evidence/SEQUENCE_CONTEXT_FALSIFICATION_01.md)

Run:
`36875022535`

Experiment result commit:
`c668ab156dbc22d02e29e5648fd0765938497d31`

Complete evidence artifact:
- ID: `11168668376`
- digest: `sha256:66bd9258f4879491d61cc45713d6b54a1f609a874fcc85e19bd64b793b67d9b4`

Decision:
`PRETRAINED_SEQUENCE_CONTEXT_SCALAR_REJECTED`

Primary quad evidence:
- TrpB partial context Spearman: -0.028251566611507472
- TrpB permutation p: 0.7202797202797203
- GB1 partial context Spearman: -0.03728501478447224
- GB1 permutation p: 0.6113886113886113

All integrity checks passed.
Deterministic replay max absolute difference: 0.0.
NucB consumed: false.
Phase 3 opened: false.

Meaning:
The preregistered scalar candidate-conditioned ESM-2 preference-shift score did not contain material higher-order residual information beyond the frozen low-order nuisance feature set and intrinsic single-mutation control.

Status:
REJECTED HYPOTHESIS INSTANCE

Scope of rejection:
This does NOT reject all sequence-aware priors. It rejects this exact scalar score, model revision, scoring rule, and gate.

## Pretrained Möbius Falsification 02

Detailed canonical record:
[evidence/PRETRAINED_MOBIUS_FALSIFICATION_02.md](evidence/PRETRAINED_MOBIUS_FALSIFICATION_02.md)

Run:
`36878521685`

Result commit:
`3ac80d29eb33f28250029131d94390a8ac40bb8d`

Complete evidence artifact:
- ID: `11171700565`
- digest: `sha256:a8cbdd632cd3acff8e18ef72e351dfc214cda1ea2bec589325120f90abb5b110`

Decision:
`PRETRAINED_MOBIUS_SCALAR_REJECTED`

Primary quad evidence:
- TrpB partial Möbius Spearman: 0.07588401649350297
- TrpB permutation p: 0.3096903096903097
- GB1 partial Möbius Spearman: -0.0046879662534248445
- GB1 permutation p: 0.9590409590409591
- cross-landscape sign consistency: false

Fresh-holdout contract:
- ranks 192..383 in every cell;
- no candidate overlap with Experiment 01 ranks 0..191;
- 192 candidates per cell;
- all integrity and hash checks passed;
- exact deterministic replay max difference: 0.0.

Meaning:
Explicit frozen ESM likelihood-based beyond-pair non-additivity did not carry transferable biological higher-order residual information after low-order nuisance control.

Status:
REJECTED HYPOTHESIS INSTANCE

Scope of rejection:
The ESM likelihood-scalar path is retired for this development track. This does not reject richer frozen embeddings/representations as a separate hypothesis class.

## Current scientific state

- Phase 1: CLOSED and preserved.
- Phase 2A/2B/2C: historical completed work.
- Phase 2D replacement gate: NOT PASSED.
- B2: frozen transfer baseline.
- Root cause: NOT YET IDENTIFIED.
- V10.0: REJECTED.
- scalar candidate-conditioned ESM-2 context score: REJECTED.
- pretrained ESM likelihood-based Möbius scalar: REJECTED.
- NucB: UNTOUCHED.
- Phase 3: CLOSED.

No architecture is currently authorized.

Any next experiment must be admitted under `SCIENTIFIC_CONSTRAINTS.md` and must answer a distinct discriminative question rather than repair a rejected score or architecture.
