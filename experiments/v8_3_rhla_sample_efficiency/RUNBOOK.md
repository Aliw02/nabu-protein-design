# NABU V8.3 RhlA Sealed Sample-Efficiency Test

This is a one-shot sample-efficiency test of the frozen NABU V8.3 architecture.

Base architecture:
- branch: `freeze/nabu-v8-3-dual-objective-router-validated`
- commit: `c8afdcd7698231d95a39ef13a3fe22b6e3f507c5`

Dataset:
- RhlA combinatorial activity landscape
- source: `sitonglab/CombinGym`
- pinned commit: `c5dccd514a9acfe1089a5f1555faea98afd78019`
- file: `Data/DMS/Clean/RhlA_clean.xlsx`
- objective: overall enzyme activity only
- selectivity and specific-activity labels are not used

Evaluation target:
- hidden triple mutants only

Training labels:
- deterministic subset of measured single, double, and triple mutants

Fixed hidden split:
- candidate mutation_count == 3
- FNV1a32(candidate_id) mod 10 in 7, 8, or 9
- hidden fitness is physically separated before Stage A

Training pool:
- mutation_count in 1, 2, or 3
- FNV1a32(candidate_id) mod 10 in 0..6

Nested experimental budgets:
- 5%
- 10%
- 20%
- 40%

The budget order is identity-only:
`FNV1a32("BUDGET|" + candidate_id)`.

All four budgets use the same frozen NABU V8.3 architecture. No tuning is allowed.

## PowerShell

From repository root:

```powershell
git fetch origin
git switch exp/nabu-v8-3-rhla-sample-efficiency-sealed
git pull --ff-only
git rev-parse HEAD

py -m pip install pandas numpy scipy openpyxl

py experiments/v8_3_rhla_sample_efficiency/00_download_rhla.py

py experiments/v8_3_rhla_sample_efficiency/00_prepare_sealed_rhla.py `
  RhlA_clean.xlsx `
  --out rhla_sample_efficiency_sealed_input
```

Stage 0 should print only schema/count/hash information. It must not print hidden triple fitness values.

Then run Stage A:

```powershell
py experiments/v8_3_rhla_sample_efficiency/01_stage_a_freeze_budgets.py `
  --input rhla_sample_efficiency_sealed_input `
  --out nabu_v8_3_rhla_sample_efficiency_stage_a
```

## Stage A qualification

Stage A runs all four nested budgets and fixes one common hidden triple-mutant pool.

The test aborts before reveal when:
- fewer than 50 hidden triple candidates are scoreable at 5%, or
- scoreable coverage is below 50% of the hidden triple test, or
- a later nested budget unexpectedly loses scoreability for any fixed candidate.

Check:

```powershell
if (Test-Path "nabu_v8_3_rhla_sample_efficiency_stage_a\STAGE_A_QUALIFICATION_ABORT.json") {
    Get-Content "nabu_v8_3_rhla_sample_efficiency_stage_a\STAGE_A_QUALIFICATION_ABORT.json"
} else {
    Get-Content "nabu_v8_3_rhla_sample_efficiency_stage_a\STAGE_A_MANIFEST.json"
}
```

If `STAGE_A_QUALIFICATION_ABORT.json` exists, do not run Stage B.

Successful Stage A freezes:
- `BUDGET_MEMBERSHIP.csv`
- `ALL_BUDGET_PREDICTIONS.csv`
- `CANDIDATE_FREEZE.csv`
- `STAGE_A_MANIFEST.json`

Only after Stage A succeeds, reveal once:

```powershell
py experiments/v8_3_rhla_sample_efficiency/02_stage_b_reveal_evaluate.py `
  --input rhla_sample_efficiency_sealed_input `
  --stage-a nabu_v8_3_rhla_sample_efficiency_stage_a `
  --out nabu_v8_3_rhla_sample_efficiency_final
```

Final outputs:
- `nabu_v8_3_rhla_sample_efficiency_final/FINAL_RESULTS.json`
- `nabu_v8_3_rhla_sample_efficiency_final/SAMPLE_EFFICIENCY_CURVE.csv`

## Frozen practical threshold

The smallest budget is considered practically sufficient relative to the frozen 40% V8.3 reference only when all are true:

1. Spearman retention >= 90%
2. Top-10 mean hidden-triple percentile is within 0.03
3. Top-50 mean hidden-triple percentile is within 0.02
4. Top-50 true Top-1% hits lose at most one hit

Each budget also reports whether V8.3 improves or protects the same-budget B3 pairwise baseline.

Negative controls:
- deterministic identity-only random-hash ranking
- 40% shuffled-label V8.3

No architecture, support, shrinkage, routing, metric, or threshold changes are allowed after hidden truth reveal.
