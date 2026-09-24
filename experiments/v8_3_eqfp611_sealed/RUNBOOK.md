# NABU V8.3 eqFP611 Sealed Higher-Order Test

This is a one-shot higher-order validation using the frozen NABU V8.3 Dual-Objective Router.

Base freeze:
- branch: `freeze/nabu-v8-3-dual-objective-router-validated`
- commit: `c8afdcd7698231d95a39ef13a3fe22b6e3f507c5`

Dataset:
- eqFP611 complete 13-site binary combinatorial fluorescence landscape
- source: `sitonglab/CombinGym`
- pinned commit: `c5dccd514a9acfe1089a5f1555faea98afd78019`
- file: `Data/DMS/Clean/eqFP611_clean.xlsx`

Candidate scope:
- mutation counts 3, 4, 5 only
- deterministic identity-only split:
  - visible: FNV1a32(candidate_id) mod 10 in 0..6
  - hidden: buckets 7..9

Stage A must build:
- at least one supported triplet entry
- at least one supported quartet entry

If either is zero, the test aborts before truth reveal and does not count as a higher-order validation.

## PowerShell

From repository root:

```powershell
git fetch origin
git switch exp/nabu-v8-3-eqfp611-sealed-higher-order
git pull --ff-only
git rev-parse HEAD

py -m pip install pandas numpy scipy openpyxl

py experiments/v8_3_eqfp611_sealed/00_download_eqfp611.py

py experiments/v8_3_eqfp611_sealed/00_prepare_sealed_eqfp611.py `
  eqFP611_clean.xlsx `
  --out eqfp611_sealed_input

py experiments/v8_3_eqfp611_sealed/01_stage_a_freeze_candidates.py `
  --input eqfp611_sealed_input `
  --out nabu_v8_3_eqfp611_stage_a
```

Stage A success requires:
- `nabu_v8_3_eqfp611_stage_a/STAGE_A_MANIFEST.json`
- no `STAGE_A_QUALIFICATION_ABORT.json`
- nonzero triplet and quartet entries in the manifest

Check:

```powershell
Get-Content nabu_v8_3_eqfp611_stage_a\STAGE_A_MANIFEST.json
```

Only after Stage A qualification passes, reveal:

```powershell
py experiments/v8_3_eqfp611_sealed/02_stage_b_reveal_evaluate.py `
  --input eqfp611_sealed_input `
  --stage-a nabu_v8_3_eqfp611_stage_a `
  --out nabu_v8_3_eqfp611_final
```

Final outputs:
- `nabu_v8_3_eqfp611_final/FINAL_RESULTS.json`
- `nabu_v8_3_eqfp611_final/REVEALED_TOP50.csv`

## Fixed pass gate

Primary endpoints:
1. hidden pooled Spearman
2. Top-10 mean hidden-fitness percentile
3. Top-50 mean hidden-fitness percentile
4. Top-50 count of true hidden Top-1% variants
5. normalized Top-1 regret

Pass requires:
- no regression against B3 on every primary endpoint
- strict improvement on at least 2 primary endpoints
- no post-reveal architecture, threshold, support, shrinkage, or routing changes
