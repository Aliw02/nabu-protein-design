# NABU V9.1 — Degree-Normalized Pair Transfer Protocol

Status: DEVELOPMENT / POST-PHASE-2 / PHASE 3 REMAINS CLOSED

Base:
- frozen V8.3 core remains unchanged: `c8afdcd7698231d95a39ef13a3fe22b6e3f507c5`
- Phase-2 first blind reveal remains immutable.
- V9.0 is preserved as a failed development mechanism.
- IRED is now a development/diagnostic dataset because its targets were already revealed.

## V9.0 failure being isolated

V9.0 factorized unseen-pair transfer increased structural coverage but degraded
ranking relative to B2 main-only.

Observed V9.0 development diagnostics:
- B2 full-test Spearman: 0.1699713866
- V9.0 pair-transfer Spearman: 0.0951920479
- B2 Top-1% hits: 3 / 42
- V9.0 Top-1% hits: 0 / 42
- B2 prediction standard deviation: about 0.98
- V9.0 prediction standard deviation: about 3.52
- V9.0 transferred-pair contributions: 20,831

The postmortem found that the absolute V9.0 correction grows strongly with the
number of transferred pair edges. The factor model expresses an edge as:

    e(i,j) = bias + u_i + u_j

Summing these edge estimates across a higher-order candidate repeats each node
propensity once for every incident edge. In a fully connected m-mutation
candidate, each u_i is repeated m-1 times.

This is an aggregation problem, distinct from the tiny cross-run numerical
sensitivity of the underdetermined least-squares parameterization.

## Single V9.1 change

Keep the V9.0 factor fit unchanged.

Keep exact frozen V8.3 pair memories unchanged.

For transferred unseen-pair edges only:

1. compute the same V9.0 transferred edge effects;
2. let E be the number of transferred edges;
3. let N be the number of mutation nodes participating in at least one
   transferred edge;
4. define mean transferred degree:

       d_bar = 2E / N

5. replace raw transferred sum S with:

       S_normalized = S / d_bar

For E = 0, the transferred contribution is zero.

For a single transferred edge, d_bar = 1, so V9.1 exactly matches V9.0.

For a complete m-node transferred graph, d_bar = m-1, which removes the
automatic m-1 repetition of each node propensity while preserving a
candidate-level contribution that grows linearly rather than quadratically.

No target-position rule, amino-acid class, dataset-name branch, learned gate,
or tuned scalar is introduced.

## Arms

1. PHASE2_V8_3_ABSTENTION
2. B2_MAIN_ONLY
3. V9_0_RAW_PAIR_TRANSFER
4. V9_1_DEGREE_NORMALIZED_TRANSFER
5. V9_1_DEGREE_NORMALIZED_PERMUTED_CONTROL

All arms use the same fit/train rows and candidate identities.

Rows with unseen main mutation identities use the fit global mean, matching
the prior V9 development contract.

## Metrics

Full 4,178-row IRED test:
- Spearman
- NDCG
- Top-1% hits / recall / enrichment
- normalized Top-1% regret

Also report:
- non-strict transfer region;
- by mutation count;
- prediction scale;
- transferred edge count;
- mean-degree normalizer distribution;
- raw versus normalized transfer magnitude;
- rank-order replay hash.

## V9.1 development decision rule

V9.1 is mechanistically promising only if all are true:

1. full-test Spearman > B2_MAIN_ONLY;
2. full-test Spearman > V9_0_RAW_PAIR_TRANSFER;
3. full-test Spearman > permuted normalized control;
4. at least one elite-selection metric improves over B2;
5. normalized Top-1% regret is not worse than B2;
6. non-strict transfer-region Spearman > B2 in that same region;
7. predictions are finite;
8. deterministic replay preserves identical rank order.

This is development evidence only.

## Boundary

A V9.1 development PASS does not reopen or repair Phase 2 and does not open
Phase 3. A new untouched external benchmark is required before a new
scientific generalization claim.
