# NABU Cross-Landscape Cross-Order Prior — Results

Date: 2026-09-28
Branch: exp/nabu-root-cause-audit
Run: 36410983124
Artifact: 10964168786
Artifact digest: sha256:a9ec41272bc8d92a5437f84cca46deba02b78415f2448efcac29da980a2cb9af

Decision:
CROSS_LANDSCAPE_CROSS_ORDER_PRIOR_UNRESOLVED

Root-cause status:
ROOT_CAUSE_NOT_YET_IDENTIFIED

NucB consumed: false
Phase 3 opened: false

## Direction 1 — GB1 triples -> TrpB quads

Residual prediction:
- Spearman: 0.9098702883
- real - shuffled: +1.3593790520
- real - source affine: -0.0364351942

Reconstructed fitness:
- transferred prior Spearman: 0.0509965080
- exact O2 Spearman: 0.0108201859
- target B2 Spearman: 0.1250627849

Predeclared direction gate: FAIL

Reasons:
- residual prediction does not beat source-affine control by >= 0.10;
- reconstructed fitness does not beat exact O2 by >= 0.10;
- reconstructed fitness is below target B2 by more than the allowed 0.02.

## Direction 2 — TrpB triples -> GB1 quads

Residual prediction:
- Spearman: 0.9612909039
- real - shuffled: +1.2888666274
- real - source affine: +0.0208380844

Reconstructed fitness:
- transferred prior Spearman: 0.2865519763
- exact O2 Spearman: 0.2027648012
- target B2 Spearman: 0.3039740784

Predeclared direction gate: FAIL

Reasons:
- residual prediction does not beat source-affine control by >= 0.10;
- reconstructed fitness improves over O2 by +0.0837871751, below the required +0.10;
- reconstructed fitness is within 0.02 of B2, but this alone is insufficient.

## Triple-to-triple replication

The previous cross-landscape triple transfer reproduced:

GB1 -> TrpB:
- residual Spearman: 0.8981811251
- shuffled control: -0.5555344949
- real - shuffled: +1.4537156200

TrpB -> GB1:
- residual Spearman: 0.8275499792
- shuffled control: 0.0000288765
- real - shuffled: +0.8275211027

## Interpretation boundary

The high cross-order residual Spearman is real relative to shuffled labels, but it is not sufficient evidence for a useful universal higher-order prior.

The source-affine cancellation control matches or exceeds most of the residual-ranking signal, and adding the transferred residual back to exact O2 does not recover target quad fitness strongly enough to pass the preregistered gate.

Therefore:
- do not promote this feature class into the NABU architecture;
- do not tune blend weights on these revealed targets;
- do not consume NucB;
- proceed only with the separately preregistered final prior-transfer ablation/stopping rule if further diagnosis is needed.
