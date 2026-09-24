# NABU V8.3 CR9114-H1 Final Phase-1 Gate

This is the final Phase-1 validation. No Phase-1 retuning is allowed after reveal.

It tests, in one fresh sealed experiment:

1. Frozen NABU V8.3 higher-order routing on hidden 4/5-mutation CR9114-H1 variants.
2. Frozen 50/50 Active Acquisition controller against passive sampling.
3. Negative controls.
4. Whether Active 20% can reach the frozen Passive 40% practical threshold.

Dataset:
- CR9114-H1 antibody binding-affinity landscape
- source: sitonglab/CombinGym
- pinned commit: c5dccd514a9acfe1089a5f1555faea98afd78019
- file: Data/DMS/Clean/bnAbs_CR9114_H1_clean.xlsx
- blob: d9895e1f7e6a25fba4c6e2aced9843f7ad4a7000

Candidate scope:
- training identities: measured 1/2/3/4/5-mutants
- hidden evaluation: 4/5-mutants only

The training-label store and hidden-test truth are physically separated from identities.

## PowerShell

From repository root:

```powershell
git fetch origin
git switch exp/nabu-v8-3-cr9114-final-phase1-gate
git pull --ff-only
git rev-parse HEAD

py -m pip install pandas numpy scipy openpyxl

py experiments/v8_3_cr9114_final_phase1/00_download_cr9114_h1.py

py experiments/v8_3_cr9114_final_phase1/00_prepare_sealed_cr9114_h1.py `
  bnAbs_CR9114_H1_clean.xlsx `
  --out cr9114_h1_final_gate_sealed_input
```

Stage 0 creates:

- `TRAINING_IDENTITIES.csv`
- `.TRAINING_LABELS.sealed.csv`
- `HIDDEN_4_5_IDS.csv`
- `.HIDDEN_4_5_TRUTH.sealed.csv`
- `SEAL_MANIFEST.json`

Do not open either sealed CSV.

Run Stage A:

```powershell
py experiments/v8_3_cr9114_final_phase1/01_stage_a_freeze_final_gate.py `
  --input cr9114_h1_final_gate_sealed_input `
  --out nabu_v8_3_cr9114_final_gate_stage_a
```

Stage A performs all of the following before hidden truth reveal:

- passive 5/10/20/40% budgets
- active 5 -> 10 -> 20 -> 40% acquisition
- identical frozen V8.3 router at every budget
- exact frozen 50/50 exploration/exploitation controller
- passive-40 shuffled-label control
- deterministic random-hash control
- common scoreable hidden 4/5-mutant pool
- Higher-Order qualification
- Top-10 / Top-50 freeze
- SHA-256 freeze

Qualification requires:

- Passive 40% triplet memory > 0
- Passive 40% quartet memory > 0
- Active 20% triplet memory > 0
- Active 20% quartet memory > 0
- at least 500 common hidden 4/5 candidates
- at least 50% hidden-pool scoreable coverage

Check:

```powershell
if (Test-Path "nabu_v8_3_cr9114_final_gate_stage_a\STAGE_A_QUALIFICATION_ABORT.json") {
    Get-Content "nabu_v8_3_cr9114_final_gate_stage_a\STAGE_A_QUALIFICATION_ABORT.json"
} else {
    Get-Content "nabu_v8_3_cr9114_final_gate_stage_a\STAGE_A_MANIFEST.json"
}
```

If the abort file exists, do not run Stage B.

If Stage A passes, reveal exactly once:

```powershell
py experiments/v8_3_cr9114_final_phase1/02_stage_b_reveal_final_gate.py `
  --input cr9114_h1_final_gate_sealed_input `
  --stage-a nabu_v8_3_cr9114_final_gate_stage_a `
  --out nabu_v8_3_cr9114_final_gate_final
```

Main outputs:

- `nabu_v8_3_cr9114_final_gate_final/FINAL_RESULTS.json`
- `nabu_v8_3_cr9114_final_gate_final/FINAL_GATE_CURVE.csv`
- `nabu_v8_3_cr9114_final_gate_final/REVEALED_TOP50.csv`

## Frozen final gates

### Core V8.3 gate at Passive 40%

Against same-budget B3:

- no regression on all five primary metrics
- strict improvement on at least 2/5 metrics

Primary metrics:

1. Spearman
2. Top-10 mean percentile
3. Top-50 mean percentile
4. Top-50 true Top-1% hits
5. normalized Top-1 regret

### Active controller gate

Active 20% is compared with Passive 40% V8.3.

Required:

- Spearman retention >= 90%
- Top-10 percentile gap <= 0.03
- Top-50 percentile gap <= 0.02
- Top-50 Top-1% hit loss <= 1

Active 10% is also measured but is not required.

## Final decision

If both the Core gate and Active-20 gate pass:

`PHASE1_FINAL_GATE_PASS`

Otherwise:

`PHASE1_FINAL_GATE_CLOSED_WITH_LIMITATION`

In either case Phase 1 ends here. Do not tune on CR9114-H1 after reveal. Move to Phase 2.
