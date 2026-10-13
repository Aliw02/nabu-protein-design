# GLP1R activity bridge: frozen result

**Decision: DO_NOT_ADVANCE**

Development-selected baseline: `RANDOM_FOREST`. Candidate: `NABU_JOINT_PAIRS`.

All models used the same lossless 44-slot sequence representation. Development used five outer edit-distance-grouped folds with three inner grouped folds; external predictions were frozen before evaluation. External primary metrics use 18 unique sequences from 19 published records.

| Method | Development RMSE | Development Spearman | External RMSE | External Spearman | External top-5 mean pEC50 |
|---|---:|---:|---:|---:|---:|
| NABU_JOINT_PAIRS | 1.499413 | 0.208251 | 0.512265 | 0.466080 | 10.334910 |
| LEGACY_B2 | 5.875818 | 0.578255 | 6.008927 | -0.004143 | 10.210880 |
| LEGACY_B3 | 629.434718 | -0.561123 | 249.365505 | 0.385293 | 10.298250 |
| JOINT_ADDITIVE | 1.684472 | 0.132973 | 0.536876 | 0.458830 | 10.334910 |
| RBF_SVR | 1.570369 | 0.211253 | 0.655959 | 0.632833 | 10.236450 |
| RANDOM_FOREST | 1.415165 | 0.311317 | 1.238530 | 0.226348 | 10.173290 |
| NEAREST_NEIGHBOR | 1.630579 | 0.189502 | 0.595893 | 0.497960 | 10.053730 |
| LENGTH_ONLY | 1.725034 | -0.529766 | 0.586536 | undefined (constant) | 9.425040 |
| TRAIN_MEAN | 1.623926 | -0.603780 | 0.503794 | undefined (constant) | 9.425040 |
| SHUFFLED_JOINT_PAIRS | 1.634287 | -0.576969 | 0.496390 | 0.356292 | 10.459870 |

Lower RMSE is better; higher Spearman and top-5 mean potency are better. These are predictions of a published cell-assay potency endpoint, not therapeutic efficacy, molecular binding validation, or new peptide experiments.

## Frozen gate

- PASS: integrity_and_full_coverage
- PASS: unchanged_deterministic_replay
- FAIL: development_rmse_gain_at_least_5pct
- FAIL: development_spearman_not_lower
- PASS: development_beats_additive
- PASS: development_beats_shuffle
- PASS: external_rmse_gain_at_least_5pct
- PASS: external_spearman_gain_at_least_005
- PASS: external_top5_mean_not_lower
- PASS: external_beats_additive
- FAIL: development_paired_mse_ci_upper_negative
- PASS: external_paired_mse_ci_upper_negative

## Paired uncertainty

Candidate-minus-champion MSE; negative values favor the candidate.

- development_oof: 0.245546; 95% CI [-0.331500, 0.603505], 32 resampling units.
- external_sequence_units: -1.271542; 95% CI [-1.651949, -0.879360], 18 resampling units.

## Boundaries

This is a standard categorical ANOVA ridge model testing a NABU restructuring hypothesis, not a novel algorithm. The legacy arms apply the original B2/B3 equations to the new adapter and are not a rerun of the complete V8.3 hierarchy. All training activity labels are allowed here, changing the question from the historical low-order-only transfer task.

Public retrospective data, one small related peptide family, incomplete chemical metadata, study-to-study assay differences, and sequence-equivalent external records limit generalization. Grouped development validation reduces near-clone leakage; it does not establish protein-family transfer. The published paper uses other protocols, so these results are not a like-for-like victory over that paper.

No post-outcome retuning is permitted. The per-case artifacts, raw source hashes, complete grids, support diagnostics, freeze records, all 19-record sensitivity metrics, and negative controls are retained. NucB and the canonical model are untouched. Competition eligibility and de novo functional validation remain separate.
