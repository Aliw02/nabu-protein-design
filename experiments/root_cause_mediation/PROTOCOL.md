# NABU Scale-vs-Context Mediation Test

Date: 2026-09-28
Branch: exp/nabu-root-cause-audit
Datasets: revealed TrpB + revealed GB1
NucB: forbidden
Phase 3: closed

## Purpose

Stage 2 supported all of:
- context-free pair non-portability (H3);
- conditional pair context dependence (H4);
- observation-scale/composition nonlinearity evidence (H5).

This test asks whether the prespecified low-order-only scale correction that transferred across both datasets is sufficient to remove H3/H4.

It does not train a new ranking architecture.

## Fixed scale candidate

robust_asinh was one of the Stage-1 preregistered transforms and was the only prespecified non-identity transform that met the Stage-2 material H5 gate on BOTH TrpB and GB1 for mutation order 3.

For each dataset its parameters are fit from low-order targets (orders 0-2) only:

z(y) = asinh((y - median_low) / MAD_low)

If MAD_low == 0, use low-order std; if that is also zero use 1.

No higher-order target may fit transform parameters.

## Correct conditional pair interaction

For a pair (a,b) and third mutation c:

base pair epistasis:
epsilon_ab = y(ab) - y(a) - y(b) + y(WT)

conditional pair epistasis on background c:
epsilon_ab|c = y(abc) - y(ac) - y(bc) + y(c)

Under the exact Mobius decomposition:
epsilon_ab|c = e2_ab + e3_abc

This is the quantity used for H4.

## Gates reused unchanged

H3 non-portability on an exact-O2-supported mutation order if either:
- O2 Spearman - O1 Spearman < 0.05
- O2 Spearman <= shuffled-O2 Spearman + 0.02

H3 portability requires BOTH:
- O2 - O1 >= 0.05
- O2 > shuffled + 0.02

H4 material context dependence on supported triples if:
- supported triples >= 100
- conditional-pair sign-flip rate >= 0.20
- e3 RMS / target SD >= 0.50

## Discriminative decisions

SCALE_CORRECTION_SUFFICIENT_FOR_H3_H4 only if on BOTH datasets under robust_asinh:
- H3 becomes portable on order 3; AND
- H4 becomes non-material.

CONTEXT_DEPENDENCE_PERSISTS_AFTER_STANDARD_SCALE_CORRECTION only if on BOTH datasets under robust_asinh:
- H3 remains non-portable on order 3; AND
- H4 remains material.

Otherwise:
SCALE_VS_CONTEXT_UNRESOLVED

Order 4 is reported but not required for the decision because H4 is defined using a third-mutation background and Stage 2's shared transform gate was met on order 3.

## Interpretation boundary

If context dependence persists after robust_asinh, this rejects the claim that this prespecified standard scale correction is sufficient. It does NOT reject every possible latent/nonlinear biological mapping.

If robust_asinh removes H3/H4, it supports scale mediation but does not by itself prove a universal observation model.

## Low-order identifiability theorem

Separately from the empirical test:

WT+single+double labels cannot uniquely identify unrestricted 3-way or higher Mobius coefficients.

For any predictor f that agrees with all observed sets |S|<=2, define:
f_lambda(S) = f(S) + lambda * h(S)
where h(S)=0 for every |S|<=2 and h(S) is nonzero for at least one |S|>=3.

Every lambda gives identical low-order observations and potentially different higher-order predictions.

Therefore higher-order prediction from low-order data always requires an additional assumption/prior or additional higher-order measurements. This theorem does NOT imply useful prediction is impossible.

## Outputs

- SCALE_CONTEXT_MATRIX.json
- TRPB_SCALE_CONTEXT.json
- GB1_SCALE_CONTEXT.json
- RUN_MANIFEST.json
- OUTPUT_HASHES.json
