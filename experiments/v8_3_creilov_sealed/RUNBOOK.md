# NABU V8.3 CreiLOV Sealed Candidate Test

This is a one-shot untouched candidate-selection test using the frozen
NABU V8.3 Dual-Objective Router.

Base freeze:
- branch: `freeze/nabu-v8-3-dual-objective-router-validated`
- commit: `c8afdcd7698231d95a39ef13a3fe22b6e3f507c5`

Dataset:
- CreiLOV combinatorial fluorescence landscape
- source: `sitonglab/CombinGym`
- pinned commit: `c5dccd514a9acfe1089a5f1555faea98afd78019`
- file: `Data/DMS/Clean/CreiLOV_clean.xlsx`

The test uses:
- visible training labels: mutation_count <= 3
- hidden candidate pool: mutation_count 4-5
- hidden fitness physically separated before Stage A
- frozen candidate budgets: Top-10 and Top-50
- deterministic shuffled-label control

## PowerShell

From repository root:

```powershell
git fetch origin
git switch exp/nabu-v8-3-creilov-sealed-candidate-test
git pull --ff-only
git rev-parse HEAD

py -m pip install pandas numpy scipy openpyxl

Invoke-WebRequest `
  -Uri "https://raw.githubusercontent.com/sitonglab/CombinGym/c5dccd514a9acfe1089a5f1555faea98afd78019/Data/DMS/Clean/CreiLOV_clean.xlsx" `
  -OutFile "CreiLOV_clean.xlsx"

py experiments/v8_3_creilov_sealed/00_prepare_sealed_creilov.py `
  CreiLOV_clean.xlsx `
  --out creilov_sealed_input

py experiments/v8_3_creilov_sealed/01_stage_a_freeze_candidates.py `
  --input creilov_sealed_input `
  --out nabu_v8_3_creilov_stage_a
```

After Stage A finishes, confirm these files exist before reveal:

```powershell
Get-Item nabu_v8_3_creilov_stage_a\STAGE_A_MANIFEST.json
Get-Item nabu_v8_3_creilov_stage_a\CANDIDATE_FREEZE.csv
Get-Item nabu_v8_3_creilov_stage_a\ALL_HIDDEN_SCORES.csv
```

Then run Stage B:

```powershell
py experiments/v8_3_creilov_sealed/02_stage_b_reveal_evaluate.py `
  --input creilov_sealed_input `
  --stage-a nabu_v8_3_creilov_stage_a `
  --out nabu_v8_3_creilov_final
```

Final outputs:

- `nabu_v8_3_creilov_final/FINAL_RESULTS.json`
- `nabu_v8_3_creilov_final/REVEALED_TOP50.csv`

## Pass gate

The test passes only when:
1. V8.3 does not regress against B3 on any primary endpoint.
2. V8.3 strictly improves at least 2 of 4 primary endpoints:
   - Top-10 mean percentile
   - Top-50 mean percentile
   - Top-50 true Top-1% hit count
   - normalized Top-1 regret

No post-reveal architecture changes, threshold changes, or reruns with modified rules are allowed for this test.
