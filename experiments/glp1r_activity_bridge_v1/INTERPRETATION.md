# Interpretation of the frozen GLP1R bridge result

**Decision: DO_NOT_ADVANCE.** The supervised joint-pair candidate did not establish a consistent advantage over strong baselines. This branch preserves a completed, rejected architecture experiment. It is not an active replacement architecture and will not receive post-outcome score repairs.

## What was actually changed and measured

The architectural intervention replaced marginal mutation-memory summation with joint regularized estimation of site effects and pair interactions. This is standard categorical ANOVA kernel ridge regression; the experiment tests its suitability as a NABU restructuring, not algorithmic novelty. A lossless 44-slot adapter was shared by every arm, so its effect is not confounded with the architecture comparison.

All 125 development GLP1R activity labels were allowed. Five outer validation folds separated connected components defined by raw-sequence edit distance <= 2, with three inner grouped folds for tuning. This is a different supervised question from NABU's historical low-order-only transfer problem. The external set contained 19 published records with 18 distinct primary sequences. All ten methods predicted every case. Prediction files, support diagnostics, champion selection, and all inner tuning scores matched exactly on a second unchanged CPU run.

## Main results

Lower RMSE is better. Higher Spearman is better. The external columns give equal weight to the 18 unique primary sequences; the complete 19-record sensitivity analysis is retained in `results/metrics.json`.

| Method | Development RMSE | Development Spearman | External RMSE | External Spearman |
|---|---:|---:|---:|---:|
| NABU joint pairs | 1.499413 | 0.208251 | 0.512265 | 0.466080 |
| Joint additive ablation | 1.684472 | 0.132973 | 0.536876 | 0.458830 |
| Random forest | **1.415165** | 0.311317 | 1.238530 | 0.226348 |
| RBF SVR | 1.570369 | 0.211253 | 0.655959 | **0.632833** |
| Nearest neighbor | 1.630579 | 0.189502 | 0.595893 | 0.497960 |
| Training mean | 1.623926 | -0.603780 | 0.503794 | Undefined |
| Shuffled-label joint pairs | 1.634287 | -0.576969 | **0.496390** | 0.356292 |

The development-selected baseline champion was random forest. The candidate's development RMSE was 5.953% worse, and its development rank correlation was lower. The paired development MSE difference was +0.245546, with a grouped-bootstrap 95% interval [-0.331500, +0.603505]. Three frozen gate criteria therefore failed.

Externally, the candidate's RMSE was 58.639% lower than the development-selected forest, and its paired MSE interval excluded zero. This is a genuine result under the specified comparison, but does not rescue the failed development gate. The training-mean and shuffled-label controls had slightly lower external RMSE than the candidate, while SVR had higher external rank correlation. The candidate was therefore not uniformly best even on the external set.

The pair arm improved on its additive ablation in both partitions. This is descriptive evidence for that particular interaction model versus that particular additive fit, not proof of mechanistic epistasis, causal residue interactions, or a robust best predictor. No separate significance claim is made for the ablation difference.

## Why a good external RMSE alone is inadequate

Post-outcome descriptive inspection found external pEC50 SD 0.487575, compared with development SD 1.479317. The external mean was 9.983042, close to the training mean of 9.856240. A constant training-mean prediction therefore achieves external RMSE 0.503794 without ranking any peptides. The shuffled arm selected alpha 100 and predicted a narrow range, 9.826246 to 9.903306. These facts explain why near-constant controls can look good on this external error metric. They were inspected only to interpret the frozen result and were not used to retune any model or gate.

The negative pooled development correlation of TRAIN_MEAN is not an implementation error: its prediction is constant within each fold but differs across folds because each training partition has a different mean. Exact external constancy is represented by null Spearman, not zero or an invented rank.

## Legacy memory behavior

The original B2 and B3 equations on the training-consensus adapter produced development RMSE 5.875818 and 629.434718, respectively; external RMSE was 6.008927 and 249.365505. Their outputs were preserved without clipping or recalibration. In this adaptation, summing marginal terms and pair residuals over many sequence changes produced extreme scales; for example, external B3 predictions ranged from -321.501475 to -3.158681 pEC50.

This establishes that the literal B2/B3 equations with this adapter are unsuitable for this GLP1R supervised task. It does **not** invalidate their reported historical low-order benchmark performance, identify the biological root cause of earlier failures, or constitute a rerun of the full V8.3 hierarchy. Beating this unstable adapter is an insufficient success criterion, which is why ordinary supervised baselines and the additive ablation were required.

## Concrete reusable contribution

- A provenance-checked, fully reproducible GLP1R activity benchmark retaining all source rows and both external aggregation views.
- An adapter that reconstructs complete sequences, including residues omitted by the source's 30-column alignment in 33 of 125 development records. The alignment was designed by the original authors for their task; omission here is a representation difference, not evidence their published study is invalid.
- Identity-grouped nested validation, an independent published study, ten preserved arms, complete hyperparameter traces, per-case errors, and support diagnostics.
- A measured negative result that rules out promoting this exact restructuring on these data.

These are research infrastructure and falsification evidence. They are not a new experimentally verified agonist, a therapeutic result, a demonstrated novel architecture, a like-for-like win over the source paper, or a competition submission.

## Next decision

Do not promote or retune this candidate using these revealed outcomes. Any future architecture study needs a distinct scientific hypothesis and an independently designated evaluation source. For GLP1R potency work, more independent measurements with complete chemical identity and assay metadata would help distinguish representation failure from study differences; adding another blend or tuning to the 18 revealed sequences would not supply that evidence.

The frozen core, NucB, and Phase 3 remain unchanged. No wet-lab work or competition submission was performed.
