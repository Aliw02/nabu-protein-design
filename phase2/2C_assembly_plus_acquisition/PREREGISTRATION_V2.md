# Phase 2C V2 — Full-Beam Maturity Probe Preregistration

## Motivation

V1 produced strong state-matched assembly proposals but the maturity gate examined only 64 proposals while the assembler's search beam contained 256.

After one successful assembly batch, V1 closed because the next truncated top-64 probe contained too few mature proposals.

## Single changed mechanism

V2 changes exactly one parameter:

- V1 maturity probe count: 64
- V2 maturity probe count: 256

The probe now covers the full CandidateAssembler beam width.

Everything else is frozen from V1:

- GB1 development landscape;
- biological reference;
- 2048 identity-only reservoir;
- 64-candidate bootstrap;
- 256 total measurements;
- batch size 8;
- target assembly order 3;
- beam width 256;
- main support >= 2;
- pair support >= 2;
- historical 50/50 acquisition;
- random control;
- pre-reveal hashes;
- hidden-truth perturbation test;
- primary metric;
- scientific PASS criteria.

## Scientific PASS criteria

Unchanged from V1.

The primary combined condition must satisfy all:

1. discovery AUC > acquisition_only;
2. final regret <= acquisition_only;
3. final Top-1% hit count >= acquisition_only;
4. at least one assembly batch is measured.

No metric is redefined based on V1 outcomes.

## Interpretation

GB1 remains development-only.

If V2 fails the same preregistered scientific gate, do not introduce a V3 by repeatedly tuning on GB1 before a broader design review.
