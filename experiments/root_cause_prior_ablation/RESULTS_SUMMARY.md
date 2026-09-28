# NABU Final Prior-Transfer Ablation — Frozen Result

Date: 2026-09-28
Branch: exp/nabu-root-cause-audit
Run: 36417988831
Artifact: 10967423198
Artifact digest: sha256:7c40554cc26d3dcf54c830f108dce3ef176f02a12caf3cde2fb4a762334452cf

Decision:
TRANSFERABLE_INTERACTION_SHAPE_SUPPORTED

Stopping action:
STOP_ROOT_CAUSE_DIAGNOSTICS_BUILD_ONE_DERIVED_ARCHITECTURE

Algebraic-prior path:
KEEP

NucB consumed: false
Phase 3 opened: false
Root-cause status: ROOT_CAUSE_NOT_YET_IDENTIFIED

## GB1 -> TrpB

O2_ONLY:
- target Spearman: 0.6853760547
- conditional residual Spearman: -0.0863073473

SHAPE_ONLY:
- target Spearman: 0.8478638424
- conditional residual Spearman: 0.5940142075

FULL:
- target Spearman: 0.8981811251
- conditional residual Spearman: 0.5748584581

Gate quantities:
- FULL conditional >= 0.30: PASS
- FULL - O2_ONLY conditional: +0.6611658054 >= 0.10: PASS
- SHAPE_ONLY target Spearman: 0.8478638424 >= 0.20: PASS

Direction gate: PASS

## TrpB -> GB1

O2_ONLY:
- target Spearman: 0.6695993506
- conditional residual Spearman: -0.2542962089

SHAPE_ONLY:
- target Spearman: 0.7648302811
- conditional residual Spearman: 0.5731807866

FULL:
- target Spearman: 0.8275499792
- conditional residual Spearman: 0.7158250808

Gate quantities:
- FULL conditional >= 0.30: PASS
- FULL - O2_ONLY conditional: +0.9701212898 >= 0.10: PASS
- SHAPE_ONLY target Spearman: 0.7648302811 >= 0.20: PASS

Direction gate: PASS

## Interpretation boundary

This ablation rejects the explanation that the previous cross-landscape triple transfer is only a scalar/nonlinear O2 cancellation shortcut.

Low-order interaction-shape features transfer substantial information between GB1 and TrpB beyond O2 alone.

This does NOT establish the full biological root cause, and it does NOT justify opening Phase 3 or consuming NucB.

Per the preregistered stopping rule:
- stop adding root-cause diagnostics;
- do not switch to sequence-aware pretrained context yet;
- derive ONE architecture from the supported interaction-shape prior;
- test that architecture on revealed development only;
- if it fails the frozen development gate, kill that architecture rather than patching it;
- only after a clean development pass may NucB be considered for the untouched blind test.
