# Clean Pair Unseen-Edge POC — Result

Status: COMPLETE
Decision: IDENTITY_ONLY_CLEAN_PAIR_TRANSFER_NOT_SUPPORTED
Phase 3: CLOSED
Scientific claim: false

## Execution

GitHub Actions run:
- 36226947062

Head:
- f8db5b61ea17e37cf5d01601217ac3821cd90bf3

Artifact:
- 10901272109
- sha256:0a79fac414c2770702f2099f7a5735ef895668c9d48d1bf4308ba8bd65653f6f

Pinned IRED source:
- aa45a2f85fb1af87b6b0e86397b3f8f292a061dc574be2536673b5bd490b0e74

All contract tests and the workflow completed successfully.

## Data actually used

Only the already-revealed IRED fit/train partition was used.

- fit rows: 3,746
- validation rows used: 0
- test rows used: 0
- WT rows: 1
- single mutants: 1,037
- double mutants: 2,708
- doubles with both singleton counterparts: 1,759
- structurally excluded doubles missing >=1 singleton counterpart: 949
- clean-pair coverage of doubles: 64.9557%
- clean-pair unique mutation nodes: 885

Five deterministic unseen-edge folds:
- fold 0: 318
- fold 1: 373
- fold 2: 353
- fold 3: 353
- fold 4: 362

Primary pooled population:
- NODE_COVERED unseen pair edges: 1,536
- ONE_NODE_COLD: 217
- TWO_NODE_COLD: 6

## Clean residual

Definition:

    epsilon(i,j) = y(i,j) - y(single_i) - y(single_j) + y_WT

Observed clean residual distribution:
- mean: 0.3098825568
- std: 1.3048353434
- min: -4.2111795426
- max: 6.0299325817
- median: 0.2836427358

## Pooled NODE_COVERED result

### ZERO_INTERACTION
- RMSE: 1.3177506707
- MAE: 1.0183018310

### TRAIN_MEAN
- RMSE: 1.2833046630
- MAE: 0.9902387341
- sign accuracy: 0.6002604167

### NODE_ADDITIVE_RELATION
- Spearman: 0.3653291423
- Pearson: 0.3742173775
- RMSE: 1.5120075128
- MAE: 1.1608492300
- sign accuracy: 0.6282552083

### PERMUTED_NODE_ADDITIVE_CONTROL
- Spearman: -0.0058374695
- Pearson: 0.0162531060
- RMSE: 1.9302652066
- MAE: 1.5339839420
- sign accuracy: 0.5019531250

## Decision checks

PASS:
- positive Spearman
- Spearman beats permuted control
- sign accuracy > 0.50
- NODE_ADDITIVE beats permuted Spearman in all 5 folds

FAIL:
- RMSE does not beat ZERO_INTERACTION
- RMSE does not beat TRAIN_MEAN

Therefore the preregistered decision is:

    IDENTITY_ONLY_CLEAN_PAIR_TRANSFER_NOT_SUPPORTED

## Interpretation

This result changes the diagnosis.

Cleaning main effects was valuable: unlike V9.0/V9.1 higher-order transfer,
the clean unseen-pair relation contains a reproducible ranking signal.

However, the minimum-norm node-additive relation is badly calibrated in
magnitude. It can rank pair residuals and predict their sign better than chance,
but its absolute residual estimates are too large/noisy to beat a constant
baseline on squared error.

Therefore the next step should not be another higher-order V9.x aggregation tweak.
The next diagnostic must determine whether the clean-pair signal can be
calibrated out-of-fold using only training-fold information, or whether the
node-additive representation itself must be replaced by a richer relational
representation.

This result does not modify Phase 2 and does not open Phase 3.
