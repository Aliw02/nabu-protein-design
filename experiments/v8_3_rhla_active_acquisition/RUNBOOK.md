# NABU V8.3 RhlA Active Acquisition Development

This is a retrospective controller-development experiment layered above the frozen NABU V8.3 architecture.

Important:
- NABU V8.3 itself is unchanged.
- RhlA truth was already revealed in the prior sample-efficiency experiment.
- This experiment therefore develops the acquisition controller retrospectively.
- If the controller is useful here, it must later be frozen and validated on a new untouched dataset.

Base architecture:
- `freeze/nabu-v8-3-dual-objective-router-validated`
- commit `c8afdcd7698231d95a39ef13a3fe22b6e3f507c5`

Frozen passive RhlA reference:
- `freeze/nabu-v8-3-rhla-sample-efficiency-40pct-threshold`
- commit `7ca494746bd6ec1a5de5e61247aa85cd396a5349`

## Controller

Initial state:
- exact same frozen passive 5% set = 118 measured variants.

Targets:
- 10%
- 20%
- 40%

At every acquisition round, unmeasured candidates are ranked using only:
1. Structural exploration:
   - support deficit across main effects and pairs
   - deficit = `1 / sqrt(1 + support)`
2. Exploitation:
   - current frozen V8.3 predicted score
   - unscoreable candidates get zero exploitation rank

Fusion:
- 50% exploration rank
- 50% exploitation rank
- deterministic candidate-id tie break

The identity selection list is written and SHA-256 frozen before labels for that selected batch are merged into the measured set.

## PowerShell

From the repository root:

```powershell
git fetch origin
git switch exp/nabu-v8-3-rhla-active-acquisition-dev
git pull --ff-only
git rev-parse HEAD

py -m pip install pandas numpy scipy openpyxl
```

The prior RhlA sealed input should still exist locally:

```powershell
Get-Item rhla_sample_efficiency_sealed_input\TRAINING_POOL.csv
Get-Item rhla_sample_efficiency_sealed_input\HIDDEN_TRIPLE_IDS.csv
Get-Item rhla_sample_efficiency_sealed_input\.HIDDEN_TRIPLE_TRUTH.sealed.csv
```

If the sealed input directory is missing, regenerate it using the prior fixed adapter:

```powershell
py experiments/v8_3_rhla_sample_efficiency/00_download_rhla.py

py experiments/v8_3_rhla_sample_efficiency/00_prepare_sealed_rhla.py `
  RhlA_clean.xlsx `
  --out rhla_sample_efficiency_sealed_input
```

Run Active Stage A:

```powershell
py experiments/v8_3_rhla_active_acquisition/01_stage_a_active_acquisition.py `
  --input rhla_sample_efficiency_sealed_input `
  --out nabu_v8_3_rhla_active_acquisition_stage_a
```

Stage A freezes:
- `ROUND_TO_10PCT_SELECTION.csv`
- `ROUND_TO_20PCT_SELECTION.csv`
- `ROUND_TO_40PCT_SELECTION.csv`
- `ACQUISITION_LOG.csv`
- `ACTIVE_BUDGET_PREDICTIONS.csv`
- `CANDIDATE_FREEZE.csv`
- `STAGE_A_MANIFEST.json`

Check Stage A:

```powershell
Get-Content nabu_v8_3_rhla_active_acquisition_stage_a\STAGE_A_MANIFEST.json
```

Then run Stage B:

```powershell
py experiments/v8_3_rhla_active_acquisition/02_stage_b_reveal_evaluate.py `
  --input rhla_sample_efficiency_sealed_input `
  --stage-a nabu_v8_3_rhla_active_acquisition_stage_a `
  --out nabu_v8_3_rhla_active_acquisition_final
```

Main outputs:
- `nabu_v8_3_rhla_active_acquisition_final/FINAL_RESULTS.json`
- `nabu_v8_3_rhla_active_acquisition_final/ACTIVE_ACQUISITION_CURVE.csv`

## Main question

Can Active 10% or Active 20% satisfy the exact practical threshold that previously required Passive 40%?

The threshold is unchanged:
- Spearman retention >= 90% of frozen Passive 40%
- Top-10 percentile gap <= 0.03
- Top-50 percentile gap <= 0.02
- Top-50 Top-1% hit loss <= 1

If Active 20% passes, that corresponds to a 50% measurement reduction relative to Passive 40%.

If Active 10% passes, that corresponds to a 75% measurement reduction relative to Passive 40%.
