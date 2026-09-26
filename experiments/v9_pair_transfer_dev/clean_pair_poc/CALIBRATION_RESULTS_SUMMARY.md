# Clean Pair Nested-OOF Calibration — Result

Status: COMPLETE
Decision: CLEAN_PAIR_AFFINE_CALIBRATION_SUPPORTED
Phase 3: CLOSED
Scientific claim: false

## Execution

GitHub Actions run:
- 36231948658

Workflow head:
- 87ebd20e4cc2b22fbf7c93659dff7dbd73476c57

Pinned IRED source:
- aa45a2f85fb1af87b6b0e86397b3f8f292a061dc574be2536673b5bd490b0e74

Artifact:
- 10903185450
- sha256:130f027bc978666bc8bcef1b0af70ed0c2ded1ba98ed8a5283f42c4dabf123a5

All calibration contracts passed and the workflow completed successfully.

## Data boundary

Only the already-revealed IRED fit/train partition was used.

- fit rows: 3,746
- validation rows used: 0
- test rows used: 0
- clean pair rows: 1,759
- primary NODE_COVERED unseen pair rows: 1,536

Calibration was learned by nested OOF prediction inside each outer training fold.
No outer holdout target entered pair fitting or calibration.

## Pooled metrics

### ZERO_INTERACTION
- RMSE: 1.3177506707
- MAE: 1.0183018310

### TRAIN_MEAN
- RMSE: 1.2833046630
- MAE: 0.9902387341
- sign accuracy: 0.6002604167

### RAW_NODE_ADDITIVE_RELATION
- Spearman: 0.3653291423
- Pearson: 0.3742173775
- RMSE: 1.5120075128
- MAE: 1.1608492300
- sign accuracy: 0.6282552083
- prediction std: 1.4137785091

### OOF_AFFINE_CALIBRATED_RELATION
- Spearman: 0.3623951488
- Pearson: 0.3698317909
- RMSE: 1.1922131390
- MAE: 0.9222608876
- sign accuracy: 0.6412760417
- prediction std: 0.4281094195

### PERMUTED_RELATION_CONTROL
- Spearman: -0.0058374695
- Pearson: 0.0162531060
- RMSE: 1.9302652066
- MAE: 1.5339839420
- sign accuracy: 0.5019531250

## Calibration stability

Outer-fold affine coefficients:

- fold 0: alpha=0.2141829874, beta=0.2763413474
- fold 1: alpha=0.2035293323, beta=0.3022006019
- fold 2: alpha=0.2089474103, beta=0.2965365916
- fold 3: alpha=0.2274682424, beta=0.3217957091
- fold 4: alpha=0.1909339617, beta=0.3111435844

All five beta values are positive and tightly clustered around ~0.30.

Calibrated RMSE beats TRAIN_MEAN in all 5 outer folds.

## Preregistered decision checks

PASS:
- calibrated positive Spearman
- Spearman preserved within 0.01 of raw
- calibrated RMSE beats ZERO_INTERACTION
- calibrated RMSE beats TRAIN_MEAN
- calibrated RMSE beats RAW_NODE_ADDITIVE_RELATION
- calibrated sign accuracy > 0.50
- beta positive in all 5 folds
- calibrated RMSE beats TRAIN_MEAN in at least 4/5 folds (observed 5/5)

Decision:

    CLEAN_PAIR_AFFINE_CALIBRATION_SUPPORTED

## Interpretation

The clean identity-only node-additive relation contains a real transferable
unseen-pair signal. Its previous failure was primarily magnitude overdispersion,
not absence of relational information.

Nested OOF affine calibration shrinks the raw prediction amplitude by roughly
70% while preserving rank order, reducing pooled RMSE from 1.5120 to 1.1922.

This is now sufficient to retain the calibrated clean-pair relation as a
candidate component for V9 development.

It does not validate higher-order transfer yet, does not repair Phase 2, and
does not open Phase 3.
