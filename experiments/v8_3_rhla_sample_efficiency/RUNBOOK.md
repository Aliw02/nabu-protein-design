# NABU V8.3 RhlA Sealed Sample-Efficiency Test

This test measures how much labeled experimental data the frozen NABU V8.3 architecture needs before candidate selection becomes practically useful.

Base architecture:
- branch: `freeze/nabu-v8-3-dual-objective-router-validated`
- commit: `c8afdcd7698231d95a39ef13a3fe22b6e3f507c5`

Dataset:
- RhlA combinatorial landscape
- source: `sitonglab/CombinGym`
- pinned commit: `c5dccd514a9acfe1089a5f1555faea98afd78019`
- file: `Data/DMS/Clean/RhlA_clean.xlsx`

Candidate scope:
- mutation counts 3, 4, 5

Fixed hidden test:
- FNV1a32(candidate_id) mod 10 in 7, 8, 9

Training pool:
- buckets 0..6

Nested visible budgets:
- 5%
- 10%
- 20%
- 40%

Budget selection is identity-only and nested by:
`FNV1a32("BUDGET|" + candidate_id)`.

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
  --out rhla_sealed_input

py experiments/v8_3_rhla_sample_efficiency/01_stage_a_freeze_budgets.py `
  --input rhla_sealed_input `
  --out nabu_v8_3_rhla_sample_efficiency_stage_a
```

## Stage A qualification

Stage A fixes a single hidden candidate pool using the smallest 5% model.

The test aborts before reveal if:
- fewer than 100 hidden candidates are scoreable at 5%, or
- 5% coverage is below 50% of the fixed hidden test.

If this file exists, DO NOT run Stage B:

```powershell
Get-Item nabu_v8_3_rhla_sample_efficiency_stage_a\STAGE_A_QUALIFICATION_ABORT.json
```

Otherwise inspect only the Stage-A manifest:

```powershell
Get-Content nabu_v8_3_rhla_sample_efficiency_stage_a\STAGE_A_MANIFEST.json
```

Expected frozen artifacts:
- `BUDGET_MEMBERSHIP.csv`
- `ALL_BUDGET_PREDICTIONS.csv`
- `CANDIDATE_FREEZE.csv`
- `STAGE_A_MANIFEST.json`

## Stage B reveal

Only after Stage A completes without qualification abort:

```powershell
py experiments/v8_3_rhla_sample_efficiency/02_stage_b_reveal_evaluate.py `
  --input rhla_sealed_input `
  --stage-a nabu_v8_3_rhla_sample_efficiency_stage_a `
  --out nabu_v8_3_rhla_sample_efficiency_final
```

Final outputs:
- `nabu_v8_3_rhla_sample_efficiency_final/FINAL_RESULTS.json`
- `nabu_v8_3_rhla_sample_efficiency_final/SAMPLE_EFFICIENCY_CURVE.csv`

## Preregistered practical threshold

The smallest budget is considered practically sufficient when, relative to the frozen 40% V8.3 reference, it satisfies all of:

1. Spearman retention >= 90%
2. Top-10 mean percentile is within 0.03
3. Top-50 mean percentile is within 0.02
4. Top-50 true Top-1% hits lose at most one hit

Each budget also compares V8.3 directly against the same-budget B3 pairwise baseline.

Negative controls:
- deterministic random hash ranking
- 40% shuffled-label V8.3
