# NABU V8.1 — GB1 Development Replication

This branch does not change the V8.1 architecture. It adds only a source-format adapter and a replication protocol.

## 1. Download the GB1 landscape

PowerShell:

```powershell
Invoke-WebRequest -Uri "https://raw.githubusercontent.com/rcronk/mutate/0861a3f1f6c58a9597fde689eb7f9a07b38e2638/ridge/data/gb1_wu2016_landscape.csv" -OutFile "gb1_wu2016_landscape.csv"
```

## 2. Convert four-letter genotypes to categorical position-state tokens

```powershell
py experiments/v8_1_gb1_replication_dev/prepare_gb1.py gb1_wu2016_landscape.csv --out GB1_Wu2016_NABU_V8_1.csv
```

Expected first-row format:

```text
mutant,DMS_score,source_variant
X1A:X2A:X3A:X4A,...
```

## 3. Run the exact V8.1 cross-fit engine

```powershell
py experiments/v8_1_crossfit_higher_order_dev/run_v8_1_crossfit_dev.py GB1_Wu2016_NABU_V8_1.csv --out nabu_v8_1_gb1_replication_results
```

## 4. Send or push the complete result directory

Expected outputs:

```text
nabu_v8_1_gb1_replication_results/
  V8_1_CROSSFIT_DEV_RESULTS.json
  V8_1_CROSSFIT_DIAGNOSTICS.json
  V8_1_CROSSFIT_REVEALED_PREDICTIONS.csv
  V8_1_CROSSFIT_TOP50.csv
```

Do not retune V8.1 from this result. We compare the frozen architecture directly against the PHOT development result.
