# NABU V8.3 — Frozen Phase-1 Architecture

## Status

**Phase 1 core architecture: frozen and closed.**

Canonical core:

`freeze/nabu-v8-3-dual-objective-router-validated`  
Commit: `c8afdcd7698231d95a39ef13a3fe22b6e3f507c5`

Phase-1 closure:

`freeze/nabu-phase1-final-cr9114-closed`  
Commit: `690494691c20edc3c9b7a6bba9a5156913173c58`

No Phase-1 architecture thresholds or routing rules should be retuned on the revealed Phase-1 datasets.

## 1. Problem formulation

A protein variant is represented as a set of discrete substitutions

[
x = {m_1, m_2, dots, m_K}.
]

Given measured variants

[
mathcal{D}_{visible} = {(x_i,y_i)}_{i=1}^{N},
]

NABU estimates fitness for unseen combinatorial variants and ranks candidates for experimental follow-up.

The Phase-1 system is deliberately non-parametric: it stores observed main effects and residual interaction memories rather than learning a large neural network.

## 2. Frozen core hierarchy

### B2 — additive main-effect memory

Global mean:

[
mu = rac{1}{N}sum_i y_i.
]

For each mutation (m), with support (n_m):

[
C(m) = (ar y_m-mu)rac{n_m}{n_m+lambda_{main}}.
]

Then

[
hat y_{B2}(x)=mu+sum_{min x}C(m).
]

### B3 — pairwise residual memory

Residual after B2:

[
r_i^{(2)}=y_i-hat y_{B2}(x_i).
]

For each pair (p):

[
R_2(p)=ar r_p^{(2)}rac{n_p}{n_p+lambda_{pair}}.
]

Then

[
hat y_{B3}(x)=hat y_{B2}(x)+sum_{psubseteq x}R_2(p).
]

B3 is the protected backbone of V8.3.

### B4 — cross-fitted triplet residual memory

B3 predictions used to create triplet targets are generated out of fold. The row being scored does not train the residual target used for its own higher-order memory.

Triplet residuals are shrinkage-weighted and require repeated support.

[
hat y_{B4}(x)=hat y_{B3}(x)+Delta_3(x).
]

### B5 — cross-fitted quartet residual memory

Quartet targets are learned from residuals after the strict cross-fitted B4 stage.

[
hat y_{B5}(x)=hat y_{B4}(x)+Delta_4(x).
]

The implementation uses support-aware confidence and partial subset coverage rather than assuming every higher-order tuple is observed.

## 3. Dual-Objective Router

The router uses **visible-data OOF diagnostics only**.

### Mode 1 — `GLOBAL_HIGHER_ORDER`

If visible OOF B4 Spearman exceeds B3 Spearman:

[
ho_{OOF}(B4) > ho_{OOF}(B3),
]

use B5 globally.

### Mode 2 — `RANK_PRESERVING_B3_TOP20_B5_RERANK`

If global OOF higher-order gain is not positive, freeze the B3 top-20% candidate region.

Higher-order reranking is allowed only inside that frozen region if the visible OOF elite diagnostic satisfies both:

- Top-50 mean true percentile improves;
- Top-50 Top-1% hit count does not decrease.

Region membership remains B3-defined, preventing cross-boundary distortion.

### Mode 3 — `B3_PROTECTED_NO_HIGHER_ORDER`

If neither global nor elite visible evidence supports higher-order use, preserve B3 exactly.

This is an explicit abstention/safety mechanism rather than a failure state.

## 4. Why the router exists

Phase-1 experiments established three different regimes:

- **Global higher-order benefit**: PHOT, TrpB, eqFP611, CR9114-H1.
- **Elite-only higher-order benefit**: GB1.
- **Higher-order suppression required**: PhoQ and CreiLOV's low-order training regime.

Therefore NABU does not assume higher-order epistasis is universally useful. The frozen architecture first asks whether the visible evidence justifies using it.

## 5. Final fresh 4/5-mutation gate

On CR9114-H1, before hidden reveal:

- Passive 40% triplet entries: 546
- Passive 40% quartet entries: 1,698
- Active 20% triplet entries: 545
- Active 20% quartet entries: 1,024
- Fixed hidden 4/5-mutant candidates: 1,673
- Hidden scoreability coverage: 100%

At Passive 40%:

- B3 Spearman = 0.9271
- V8.3 Spearman = 0.9387
- 4/5 primary metrics strictly improved
- 5/5 primary metrics had no regression

This passed the frozen core gate.

## 6. Optional Active Acquisition layer

The Active Acquisition controller is **outside** the V8.3 core.

Frozen development checkpoint:

`freeze/nabu-v8-3-rhla-active-acquisition-dev-2x`

Policy:

[
AcquisitionScore =
0.5 	imes ExplorationRank +
0.5 	imes ExploitationRank.
]

Exploration favors poorly supported mutation/pair structure. Exploitation favors candidates ranked highly by the current frozen V8.3 model.

The controller follows a measure-update-measure loop, but Phase-1 evidence showed **non-monotonic budget behavior**:

- RhlA retrospective development: Active 20% reached the Passive 40% practical threshold.
- CR9114-H1 fresh gate: Active 10% reached Passive 40%, but Active 20% failed the preregistered target.

Therefore the controller is preserved as an experimental outer policy, not part of the universally validated core.

## 7. Anti-leakage protocol

Final sealed experiments use:

- deterministic identity-only splits,
- physically separated identity and truth files,
- SHA-256 manifests,
- frozen candidate rankings before truth reveal,
- negative shuffled-label controls,
- no post-reveal retuning.

Historical runbooks remain unchanged to preserve auditability.

## 8. Phase-2 boundary

Phase 2 will build **above** the frozen V8.3 core.

The first goal is not another static benchmark. It is a closed-loop combinatorial design system that:

1. starts from a small measured set,
2. assembles or selects candidate variants,
3. chooses which candidates to measure,
4. receives only those results,
5. updates memory,
6. repeats under a fixed experimental budget.

Phase 2 may revise the outer design/acquisition policy, but the Phase-1 V8.3 core should remain available as the frozen reference baseline.
