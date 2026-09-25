# V9.0 numerical determinism audit

Status: LABEL-INDEPENDENT ENGINEERING FIX

V9.0 V1 produced the same rankings, decision, and scientific metrics in two
successful CI runs, but raw factorized prediction floats differed at
sub-picounit scale across runners.

Observed maximum absolute cross-run difference:
- V9_PAIR_TRANSFER: 4.476419235288631e-13
- V9_PAIR_TRANSFER_PERMUTED_CONTROL: 7.185363415374013e-13

The PHASE2_V8_3_ABSTENTION and B2_MAIN_ONLY arms were bit-identical.
The V9 and permuted-control rank orders were also identical.

This is classified as LAPACK/BLAS floating-point noise in the minimum-norm
least-squares solve, not a scientific-model difference.

V2 makes only label-independent numerical-stability changes:
- same design matrix;
- same residual targets;
- same numpy least-squares solve and rcond;
- canonicalize fitted coefficients to 12 decimal places;
- CI uses single-thread BLAS/OpenMP;
- all V9.0 arms, metrics, and decision rules remain unchanged.

Preserved V1 CI records:
- run 36195301072, artifact 10890505328
- run 36195319675, artifact 10889661516

Both V1 runs concluded:
V9_PAIR_TRANSFER_NOT_YET_SUPPORTED.

The V2 rerun is an engineering reproducibility check and does not create a new
scientific hypothesis.


## V2 cross-run outcome

V2 preserved the same scientific metrics and rank orders but did not achieve
bitwise cross-run equality.

Runs:
- run 36195648704 (coefficient canonicalization)
- run 36195665646 (single-thread BLAS)
- run 36195668525 (same stabilized head plus audit commit)

Across stabilized runs the V9 rank order remained identical, but raw prediction
differences remained on the order of 1e-12 to 1e-11.

The factor system has 1299 parameters (bias + 1298 mutation nodes) and design
rank 1290. Therefore the node-factor representation is underdetermined and
numerically sensitive in its minimum-norm parameterization.

This is now treated as a property of the V9.0 factorization design rather than
a CI-only defect. No further numerical patching will be used to rescue V9.0.
