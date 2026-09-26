# Clean Pair Calibration POC

Status: POST-PHASE-2 DEVELOPMENT DIAGNOSTIC
Phase 3: CLOSED
Purpose: perform one final calibration test before deciding whether to retain or abandon the identity-only node-additive pair representation.

## Starting evidence

The clean-pair unseen-edge POC established that, on 1,536 NODE_COVERED unseen double-mutant edges:

- NODE_ADDITIVE_RELATION Spearman = 0.3653291423
- Pearson = 0.3742173775
- sign accuracy = 0.6282552083
- RMSE = 1.5120075128

Controls:
- ZERO_INTERACTION RMSE = 1.3177506707
- TRAIN_MEAN RMSE = 1.2833046630
- PERMUTED relation Spearman ≈ -0.00584

Therefore the representation contains reproducible rank/sign signal, but its magnitude is miscalibrated.

## Question

Can a calibration learned only from training-fold information preserve the clean-pair ranking signal while reducing absolute error below the constant baselines?

## Frozen data semantics

Use only IRED:

    set == train AND validation != True

Validation and test partitions remain unused.

Clean residual remains:

    epsilon(i,j) = y(i,j) - y(single_i) - y(single_j) + y_WT

Use the same outer five unseen-edge folds from the previous POC.

Primary evaluation population remains:

    NODE_COVERED unseen double-mutant edges

## Single calibration mechanism

For each outer fold:

1. Outer-train = four folds.
2. Outer-holdout = one fold.
3. Fit the raw NODE_ADDITIVE_RELATION only on Outer-train.
4. To learn calibration without leakage, generate raw predictions for Outer-train by an inner five-fold unseen-edge cross-fit performed only inside Outer-train.
5. Keep only inner predictions whose two mutation nodes are represented in that inner training split.
6. Fit one ordinary least-squares affine map:

       epsilon = alpha + beta * raw_prediction

   using those inner OOF predictions and their true clean residuals.
7. Apply the frozen alpha/beta to raw NODE_ADDITIVE_RELATION predictions on the outer holdout.

No regularization, clipping, grid search, threshold search, isotonic fitting, or post-hoc scalar tuning is allowed.

## Arms

A. ZERO_INTERACTION
B. TRAIN_MEAN
C. RAW_NODE_ADDITIVE_RELATION
D. OOF_AFFINE_CALIBRATED_RELATION
E. PERMUTED_RELATION_CONTROL

## Metrics

Primary pooled NODE_COVERED outer holdout edges:
- Spearman
- Pearson
- RMSE
- MAE
- sign accuracy

Also report:
- alpha and beta by outer fold;
- number of inner OOF rows used for calibration;
- raw vs calibrated prediction standard deviation;
- fold-level metrics;
- deterministic prediction hashes.

## Decision rule

Calibration is considered sufficient to retain the identity-only pair representation only if ALL are true:

1. calibrated Spearman > 0;
2. calibrated Spearman is not lower than raw Spearman by more than 0.01;
3. calibrated RMSE < ZERO_INTERACTION RMSE;
4. calibrated RMSE < TRAIN_MEAN RMSE;
5. calibrated RMSE < RAW_NODE_ADDITIVE_RELATION RMSE;
6. calibrated sign accuracy > 0.50;
7. beta > 0 in all five outer folds;
8. calibrated RMSE beats TRAIN_MEAN in at least 4 of 5 outer folds.

If this fails, stop node-additive calibration work. The next architecture must replace the identity-only pair relation with a richer relational representation.

## Boundary

A PASS retains this representation as a candidate component for V9 development only.
A FAIL closes this representation family.
Neither outcome changes Phase 2 or opens Phase 3.
