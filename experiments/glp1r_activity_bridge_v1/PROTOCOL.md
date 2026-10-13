# NABU GLP1R Activity Bridge V1

Protocol drafted: 2026-10-12; finalized before execution on 2026-10-13. This document and executable configuration must be committed before activity-data model fitting or comparative outcome calculation. Synthetic implementation-integrity fits are permitted before the commit.

## Authorization and scope

The user explicitly requested a measurable architectural restructuring for a useful contribution to GLP1R work. This is a bounded, supervised **peptide activity prediction** experiment. It changes the scientific question and allowed inputs for this experiment only. It does not claim to solve the previous target-WT/singles/doubles-only higher-order transfer problem. It does not reopen Phase 3, consume NucB, rewrite historical gates, or authorize changes to the frozen core. The existing frozen-embedding branch is not modified or run. This branch contains the new work only.

## Scientific question and competing hypotheses

Can joint, regularized estimation of residue-state main effects and pair interactions improve measured GLP1R potency prediction over NABU's frozen marginal memory equations and strong same-information supervised baselines?

- H1: joint pair estimation supplies useful information beyond joint additive estimation, including on the external study.
- H0a: the main benefit, if any, is ordinary additive regression rather than interactions. Reject H1 if the pair arm cannot improve on the additive ablation.
- H0b: apparent success is near-duplicate interpolation or length recognition. Use edit-distance grouped validation, length-only and nearest-neighbor controls.
- H0c: high apparent performance can occur without sequence/activity association. Fit the same pair pipeline after one deterministic training-label permutation.
- H0d: any development gain fails to transfer across assays. Freeze external predictions before revealing external targets to the evaluator.

## Sources and dataset roles

Puszkarska et al., Nature Chemistry (2024), DOI 10.1038/s41557-024-01532-x. Official authors' source: https://github.com/amp91/PeptideModels at commit `0a23c0bb86ed251acb671db09295aeac7c9c11cf`. The upstream Apache-2.0 license is retained under `third_party/`.

- `training_data.xlsx`: 125 measured sequences. Used for nested grouped development validation and final fitting.
- `reference_data.xlsx`: 19 Day et al. records. Used once for retrospective external evaluation after prediction freeze. Two records share a primary sequence, giving 18 distinct sequences. Do not claim prospective wet-lab validation or globally unseen public data.
- The source paper and public aggregate benchmark results were read before this protocol. Individual external targets have not been used for feature choice, fitting, hyperparameter selection, or method comparison. Data-loading/schema checks are not an outcome blind; public targets are not cryptographically inaccessible to the operator.
- No synthetic activity labels, docking scores, or pretrained-model predictions are treated as experimental truth.

## Information contract

Allowed: training peptide identities, their published alignment, all 125 training GLP1R potency values regardless of mutation order, and external peptide identities. GCGR labels are not model inputs. External GLP1R labels may enter only the evaluator after frozen predictions. No target labels from old NABU benchmarks are loaded.

Target is `-EC50_LOG_T2` for development and `-ec_50_log_GLP1R` externally. Higher pEC50 means stronger reported potency. Preserve supplied rounded log values. Verify raw-pM/log-molar consistency within 0.01 on training records. Do not impute, clip, invent dose-response maxima, or interpret reported EC50 as therapeutic efficacy. The source does not provide sufficient chemistry metadata or Emax for those claims.

## Lossless sequence adapter

The published training alignment is 30 columns and removes residues from some complete sequences. For each record, locate the ungapped published alignment as a contiguous substring in the full sequence. Preserve its N-terminal prefix and C-terminal suffix in separate slots. Require exact reconstruction of every raw sequence; fail if an ambiguous/missing mapping cannot be represented. Prefix slots are right aligned, suffix slots left aligned. Use the maximum prefix/suffix lengths in development identities, plus the 30 core columns; external overflow is a preparation failure, not silently cropped data. This identity-only schema is shared by all methods and frozen before fitting. This is a declared label-free, transductive use of development sequence identities for column layout, not test-target fitting.

All models see the same lossless categorical slots. Gaps are categorical states. No candidate is removed because a mutation or pair is unsupported. Report legacy memory support and unseen-state rates. Do not encode artificial chemical modifications or treat same-sequence records as proven experimental replicates.

## Architectural intervention and novelty boundary

The candidate replaces marginal lookup/summed residual composition with a **jointly fitted additive-plus-pair categorical ANOVA kernel ridge model**. This is a standard statistical formulation used to test a NABU restructuring hypothesis, not a claim of a newly invented algorithm.

At each fit, center each position's one-hot states using training frequencies. For samples x,z, let k_i be the inner product at position i. Main kernel is sum_i k_i. Pair kernel is ((sum_i k_i)^2 - sum_i k_i^2)/2. Normalize each nonzero kernel by its training mean diagonal. Candidate kernel is their sum with equal fixed prior scale. Center the final train/cross kernels for an unpenalized intercept. Select only ridge regularization by inner grouped validation. The additive ablation removes the pair kernel and otherwise uses the identical fitting procedure. There are no post-result blend coefficients, biological hotspot choices, or external-label weights.

## Baselines and controls

- LEGACY_B2 and LEGACY_B3 import `fit_base` / `score_base` from the unchanged canonical V8.1 hierarchy module used by V8.3. The adapter defines a training-only consensus reference and string mutation tokens including gap states. This is an application of frozen equations to a new input domain, **not a rerun of the entire V8.3 triplet/quartet router**.
- JOINT_ADDITIVE: normalized, centered single-position kernel ridge, the pair ablation.
- RBF_SVR: one-hot lossless features, fixed epsilon 0.1, C and gamma grid in JSON; gamma is divided by the slot count.
- RANDOM_FOREST: 200 trees; training-only grid for leaf size and feature fraction.
- NEAREST_NEIGHBOR: average labels of equally nearest lossless categorical sequences.
- LENGTH_ONLY: ridge on core nongap count, N-prefix length, and C-suffix length, with training-only scaling.
- TRAIN_MEAN.
- SHUFFLED_JOINT_PAIRS: same candidate and tuning with shuffled training labels. This is one falsification control, not a full permutation p-value.

Each hyperparameter is chosen using only the current fit partition. No pretrained paper model is used on its own training data.

## Splits and fair evaluation

Connect development sequences whose raw-sequence Levenshtein distance is at most 2. Connected components define groups, so no near-clone pair crosses a split. Use deterministic 5-fold GroupKFold outside, 3-fold GroupKFold inside. Each of 125 development cases receives exactly one out-of-fold prediction per arm. Final external fitting uses all 125 rows; hyperparameters are selected by 3-fold grouped CV on those rows. Use the same splits for every arm. Freeze group membership and fold manifests before fitting.

External predictions retain all 19 records. Primary external metrics average targets and predictions over identical primary sequences (18 units), giving equal weight to sequence-identical records; raw 19-record metrics and their target disagreement are a mandatory sensitivity report. Unknown chemical differences are a limitation, not silently resolved by averaging. External sequence identities must not duplicate development identities.

## Metrics and gate

Primary: RMSE in pEC50 units. Secondary: MAE, Spearman, top-5 mean measured pEC50, top-5 recall of the true best five, full coverage, per-case error, and runtime. Also report grouped-fold metrics and activity ranking within the better-than-1-nM subset as descriptive only when at least five cases are present; this subset never controls training or candidate removal.

Select the baseline champion by smallest pooled development OOF RMSE from the preregistered baseline list; break ties lexically. This choice is frozen before external labels enter the evaluator.

ADVANCE_TO_INDEPENDENT_REPLICATION requires all:

1. Full coverage, exact identity split checks and freeze/hash checks pass.
2. Candidate development RMSE improves at least 5% over the champion, development Spearman is not lower, and candidate RMSE is below the shuffled control and additive ablation.
3. External candidate RMSE improves at least 5% over the development-selected champion, external Spearman improves by at least 0.05, and external top-5 mean potency does not decrease.
4. External candidate RMSE is below the additive ablation.
5. The upper 95% paired-bootstrap bound on candidate-minus-champion MSE is negative, in both development and external evaluation. Development bootstrap resamples sequence components. External bootstrap resamples the 18 unique sequence units; relatedness within this one family limits generalizability.

If any efficacy criterion fails: DO_NOT_ADVANCE. If integrity fails: INVALID, with no efficacy interpretation. Even ADVANCE means only replicate this supervised predictor on an independent dataset; it does not authorize competition sequences or establish de novo agonism.

## Kill rule, determinism and evidence

No retuning or new score weights after outcomes. Preserve all failed arms. Reproduction is allowed only with unchanged protocol and code; fixes to verified implementation defects must be logged, and any outcome-exposed rerun is development, not fresh external evidence. Use the JSON seed, deterministic splits, single-worker tree fitting and fixed BLAS thread counts. Run implementation-integrity tests before fitting. Freeze the full predictions and selections before evaluation. Save raw source hashes, normalized data, split files, tuning traces, timings, software versions, prediction hashes, per-case outcomes, metrics and DECISION.json.

No ESM prior is added and none of the rejected ESM scalar paths is repaired. No NucB file is accessed. Competition submission remains a separate de novo and functional-validation question.

## Execution details fixed before outcomes

Use Python 3.12 and the pinned requirements. When an inner-CV score ties, choose the first entry in the JSON grid; the baseline-champion tie breaker remains lexical. Use seed `20261012 + outer_fold` for each development fit and `20261012 + 1000` for each final fit. The shuffled arm permutes the fit partition's labels once before its inner CV. Tree fits use that same per-partition seed for each grid choice. No labels are standardized, clipped, or filtered.

Top-five ties are broken by candidate identifier; sequence-aggregated external identifiers are the sorted source identifiers joined with `|`. Constant predictions have undefined (null) Spearman; they cannot satisfy a rank gate. The better-than-1-nM descriptive subset means pEC50 > 9. Grouped-bootstrap samples draw the same number of components with replacement and include all rows in each drawn component; the per-resample MSE difference is averaged over those rows. The external bootstrap uses one unit per unique primary sequence. Percentile 2.5% and 97.5% bounds define the 95% interval.

Before opening external outcomes, perform one unchanged deterministic replay and require identical development and external prediction artifacts. Preserve that check and commit the full prediction freeze. Evaluation creates an INVALID decision on an integrity exception. No competing architecture development is started in parallel; the previous unresolved frozen-embedding result is preserved without rerunning or modifying it.
