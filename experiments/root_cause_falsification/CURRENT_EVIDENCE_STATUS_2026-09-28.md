# NABU Root-Cause Falsification — Current Evidence Status

Date: 2026-09-28
Branch: exp/nabu-root-cause-audit

## Global status

ROOT_CAUSE_NOT_YET_IDENTIFIED

Phase 3: CLOSED
NucB consumed: false

## Canonical runs

Stage 1:
- run 36398567359
- artifact 10958833001
- digest sha256:5f5ce486010a73d54807ad8a6f343a1b8ace19d233526e89c674aec9ab82b1ff

Stage 2:
- run 36399224194
- artifact 10960356525
- digest sha256:56b69ceac6dc27c1286b5854745eadcc11df9911a63979b1ce39d4c9c4abd06e

Scale-vs-context mediation:
- run 36401928493
- artifact 10960542278
- digest sha256:aec5080ed58720a9380fab3e000be882ca7559ed2ef8a846dbee03edc2fac9e1

Stage 3:
- canonical successful run 36402076016
- artifact 10961145513
- digest sha256:0724e9fe003336cba9afa7ce2682bfcdd84a3a8bfa0b8218b84eebad761c29eb

## What is actually rejected

### Implementation / data semantics as the main explanation
REJECTED by Stage 1.

Exact O1/O2 reconstruction invariants, mutation parsing, support bookkeeping, metric self-tests, split counts, and source hashes passed on TrpB and IRED.

### Basic low-order representation-direction failure
REJECTED as the primary explanation by Stage 1.

TrpB:
- exact singleton vs B2 main Spearman: 0.7998737410
- exact pair e2 vs B3 pair contribution Spearman: 0.7093196656

IRED:
- exact singleton vs B2 main Spearman: 0.6707548
- exact pair e2 vs B3 pair contribution Spearman: 0.7251674

The learned effects are strongly shrunk in scale, but their directions are materially recovered.

## What is supported across two revealed landscapes

### Context-free pair non-portability
SUPPORTED on TrpB and GB1.

Order 3:
- TrpB exact O2 - O1 Spearman: -0.8800378347
- TrpB exact O2 - shuffled O2: -0.6914913265
- GB1 exact O2 - O1: -0.3376840878
- GB1 exact O2 - shuffled O2: -0.1585081425

Order 4 also satisfies the preregistered non-portability criterion on both datasets.

Interpretation:
an exact pair effect measured in a double cannot be treated as a portable context-free additive contribution at higher mutation order on these two landscapes.

### Conditional pair context dependence
SUPPORTED on TrpB and GB1.

The exact conditional pair interaction is:
epsilon_ab|c = y(abc) - y(ac) - y(bc) + y(c)
             = e2_ab + e3_abc

Observed scale:
- TrpB sign-flip rate: 0.3630964433
- TrpB e3 RMS / target SD: 5.9524430494
- GB1 sign-flip rate: 0.2695160109
- GB1 e3 RMS / target SD: 1.6838263644

## Scale nonlinearity: contributory, not sufficient

robust_asinh was fit from low-order targets only.

It materially improves exact O2 ranking on order-3 rows:
- TrpB Spearman gain: +0.3091844113
- GB1 Spearman gain: +0.4362229667

But it does NOT remove pair non-portability or conditional context dependence.

After robust_asinh:
TrpB:
- sign-flip rate: 0.3859080786
- e3 RMS / target SD: 2.3076629841
- H3 non-portability remains true
- H4 material remains true

GB1:
- sign-flip rate: 0.3312296119
- e3 RMS / target SD: 1.0594111195
- H3 non-portability remains true
- H4 material remains true

Therefore the tested standard observation-scale correction is not a sufficient explanation.

This does not reject every possible latent/nonlinear biological mapping.

## Simple cancellation law is not universal

A scalar affine cancellation law fit on supported triples and frozen before applying to quads:

TrpB:
- quad Spearman improvement over raw O2: +0.1005185260
- passes the preregistered transfer gate

GB1:
- quad Spearman improvement over raw O2: 0.0
- fails the transfer gate

Therefore one global scalar cancellation law is not supported across landscapes.

## Material fourth-order residual remains after scale correction

After robust_asinh and measured O3 reconstruction:

- TrpB e4 RMS / target SD: 2.8157971835
- GB1 e4 RMS / target SD: 1.6952732515

Both are above the preregistered material threshold 0.50.

The Stage-3 observation-scale promotion gate therefore fails.

## Coverage bias

UNRESOLVED / not a consistent root explanation.

TrpB support effects were below the preregistered 0.10 material threshold, while GB1 showed larger support effects and order-4 direction did not match TrpB.

## Important correction to Stage-2 H6

Stage 2 reported H6_LOW_ORDER_INFORMATION_INSUFFICIENCY = REJECTED because a random-forest model predicted higher-order residuals with high cross-validated Spearman.

That result MUST NOT be interpreted as proof that low-order labels alone are sufficient.

Reason:
the RF features were low-order-derived, but the RF was trained on higher-order residual TARGETS in its training folds.

It proves:
- higher-order residuals are strongly structured and learnable after observing higher-order labels.

It does NOT prove:
- that the residual law can be learned from WT/single/double labels alone.

For zero-shot low-order sufficiency, Stage-2 H6 is therefore non-adjudicating.

## Mathematical identifiability constraint

Without an additional structural prior, WT + singles + doubles cannot uniquely determine unrestricted third-order or higher Mobius coefficients.

For any function f matching all observed sets |S| <= 2, define:
f_lambda(S) = f(S) + lambda * h(S)
where h(S)=0 on every |S|<=2 but h(S) is nonzero on at least one |S|>=3.

Every lambda produces identical low-order evidence and different higher-order behavior.

This is a mathematical constraint, not a dataset-specific model failure.

It does NOT imply useful zero-shot prediction is impossible.
It means useful zero-shot prediction requires a justified inductive prior or additional higher-order measurements.

## Current engineering implication

Do NOT add:
- another scalar pair correction,
- another order-specific additive e3/e4 memory,
- another global cancellation coefficient,
- a local reranker presented as the root fix,
- a transform-only fix.

The next discriminative question is whether an external/pretrained sequence-context prior can predict the context-dependent interaction change using ONLY low-order supervised labels.

Until that is tested cross-dataset, no replacement architecture is frozen and NucB stays untouched.
