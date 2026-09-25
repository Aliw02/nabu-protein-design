# NABU Phase 2D — Frozen Blind Execution Preregistration

**Version:** NABU_PHASE2D_BLIND_V1  
**Date:** 2026-09-25  
**Depends on:** NABU_PHASE2C_MULTILANDSCAPE_V3 (FROZEN PASS)

## Scientific rule

Phase 2D is evaluation only. No Phase-2C architecture, support threshold, router rule,
assembly rule, batch rule, acquisition rule, or Phase-1 V8.3 mathematics may be
retuned after blind-data contact.

All failed runs are preserved.

## 2D-A — FLIP AAV sampled/random active optimization

Published protocol source:
Greenman, Amini & Yang, PLOS Computational Biology (2025),
"Benchmarking uncertainty quantification for protein engineering".

Pinned released-code reference:
- repository: microsoft/protein-uq
- commit: 5e7b2b9cd219805eaabe73b21d4d9955cf882448
- AAV source archive:
  https://raw.githubusercontent.com/microsoft/protein-uq/5e7b2b9cd219805eaabe73b21d4d9955cf882448/data/aav/splits.zip

Published campaign contract:
- AAV sampled/random training pool only;
- initial observed subset = 10% of that training pool;
- 3 folds / initialization seeds: 0, 1, 2;
- 5 total log-spaced training-set sizes from 10% through 100%;
- identical initial subsets/checkpoints for every compared method;
- primary metric = AUC of global Top-100 recovery vs measurements/fraction observed;
- report mean, standard deviation, per-fold values, and paired deltas.

Competitors:
- frozen NABU Phase-2 fixed-pool acquisition path;
- greedy;
- UCB;
- Thompson Sampling;
- random.

### Fairness boundary for assembly

The published AAV BO benchmark is a fixed-candidate-pool comparison.
No method may assay a sequence outside that common pool.

Therefore Phase-2C generative assembly is not allowed to enlarge the AAV pool.
2D-A evaluates frozen NABU ranking/acquisition under fixed-pool fairness.
Assembly evidence remains the preregistered 2C V3 development result and is
reported separately in the final Phase-2 verdict.

### AAV representation capability gate

Before target values are used, identity-only preflight must report:
- train/test counts;
- sequence lengths;
- character alphabet;
- duplicate identities;
- fraction representable by the frozen NABU substitution/deletion token grammar.

No target column may be read by the preflight code.

If the complete published AAV training pool cannot be represented without
changing NABU input semantics, record ADAPTER_SCOPE_FAIL. Do not silently
discard unsupported identities and do not redefine the benchmark after seeing
fitness outcomes.

## 2D-B — FLIP2 IRED two-to-many

Source:
- official FLIP2 download:
  https://flip.protein.properties/assets/splits/ired/two_to_many.csv.gz
- first-touch SHA256 is recorded and becomes mandatory for all reruns.

Frozen split:
- fit/train: official `set == train` rows with `validation != True`
  (3746 rows in the first-touch official file);
- validation: official `validation == True` rows are held out from fitting
  (662 rows) and are not used to tune Phase-2D;
- test: official `set == test` higher-order variants (4178 rows);
- official set/validation assignments are authoritative.

Primary metrics:
- Spearman rank correlation;
- NDCG.

Secondary metrics:
- Top-k enrichment;
- normalized regret;
- Top-1% recall where defined.

Published FLIP2 comparison values are contextual baselines; direct comparisons
are made only where metric definitions are identical.

### IRED representation capability gate

Identity-only preflight must verify fixed-length standard-amino-acid sequences
and derive mutation identities relative to the unique/reference sequence using
only sequence/set metadata. If reference derivation is ambiguous, execution
halts before target use and records ADAPTER_SCOPE_FAIL.

## Source-contact rule

The first blind-data operation is identity-only preflight.
It may download the official files but must not access the target column.
It records complete-file SHA256 values, source URLs, archive members, schema,
identity counts, sequence-length distributions and alphabets.

After first-touch hashes are recorded, no different source bytes may replace
them in the same blind version.

## Post-reveal bug rule

Only demonstrably label-independent implementation bugs may be fixed.
Every fix requires:
1. preserving the failed run;
2. an AUDIT_NOTE;
3. a new version;
4. a complete rerun.

Scientific retuning is prohibited.

## IRED full-test scoreability rule

Before any target values are read, the identity-only preflight must compute
frozen-V8.3 structural scoreability on every official test identity.

A test identity is structurally scoreable iff:
- every mutation token exists in the fit/train main-effect identity support; and
- at least one internal mutation pair exists in the fit/train pair identity support.

No test row may be silently dropped.

If structural scoreability is below 100%, a full-test handling rule must be
versioned and frozen before target reveal. Target values may not be used to
choose that rule.

## IRED secondary-metric definitions

Let `N` be the complete official test set and
`K = ceil(0.01 * N)`.

- **Top-1% recall:** intersection size between predicted top-K and true top-K,
  divided by K.
- **Top-1% enrichment:** observed true-top-K fraction among predicted top-K
  divided by the random expectation K/N. Equivalently
  `hits * N / K^2`.
- **Normalized regret at top-1% budget:** let `b` be the best true target among
  predicted top-K, and let `y_min/y_max` be the global test extrema.
  Report `(y_max - b) / (y_max - y_min)`; 0 is optimal.
- **NDCG:** match FLIP2 baseline code: shift test targets by their minimum to
  make relevance nonnegative, then use scikit-learn `ndcg_score` with the
  frozen NABU prediction as the ranking score.
- **Spearman:** scipy `spearmanr(target, prediction)`.

All metrics use the complete official test set. No scoreable-only metric may be
substituted for the full-test primary result.


## IRED frozen full-test evidence backoff

The identity-only scoreability gate found that strict frozen-V8.3 scoreability
does not cover the complete official test domain. Before any target reveal,
the following `IRED_EVIDENCE_BACKOFF_V1` policy is frozen:

- all main effects supported + at least one supported internal pair:
  use unchanged `V8_3_ADAPTIVE_ROUTER`;
- all main effects supported + no supported internal pair:
  use unchanged `B2_ADDITIVE`;
- one or more unseen main mutation identities:
  use the fit/train global mean as an abstention score.

Equal scores retain stable source order. The full 4,178-row test set remains the
primary metric domain. The strict-scoreable subset is diagnostic only.

This changes no learned parameter, memory, threshold, router rule, or Phase-1
mathematics and is frozen from identity evidence only.

## IRED V2 full-test abstention policy — frozen before target reveal

Identity-only scoreability preflight on the first-touch official source found:

- official test rows: 4178;
- frozen-V8.3 structurally scoreable: 131 (3.135471517%);
- structurally unscoreable: 4047;
- rows with unseen main identity: 530;
- rows without any supported pair identity: 4045.

No target values were read to obtain these counts.

To preserve the frozen V8.3 scoreability semantics while still evaluating the
complete official test set:

1. score structurally scoreable rows using the unchanged
   `V8_3_ADAPTIVE_ROUTER`;
2. do not invoke B2/B3/B4/B5 as a new fallback for rows that frozen V8.3 marks
   unscoreable;
3. assign every unscoreable row the same deterministic abstention score,
   strictly below the minimum finite router score among scoreable test rows;
4. compute primary Spearman and NDCG on all 4178 rows;
5. report structural scoreability/abstention coverage beside every result;
6. scoreable-only metrics may be reported only as diagnostics and may never
   replace the full-test primary metrics.

The abstention floor affects rank only. Its numeric value is
`min(scoreable_router_score) - 1.0`.

This policy is an evaluation wrapper for frozen abstention, not a scientific
retune of NABU.
