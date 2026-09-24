# V8 Adaptive Higher-Order — Development Run

This branch is intentionally not frozen. It is a development experiment on PHOT after V7.

Run from the repository root:

```powershell
py -m pip install pandas numpy scipy
py experiments/v8_adaptive_higher_order_dev/run_v8_dev.py PHOT_CHLRE_Chen_2023.csv --out nabu_v8_dev_results
```

The result should produce:

```text
nabu_v8_dev_results/
  V8_DEV_RESULTS.json
  V8_DEV_REVEALED_PREDICTIONS.csv
  V8_DEV_TOP50.csv
  V8_DEV_MEMORY_DIAGNOSTICS.json
```

Send the complete terminal output or the whole `nabu_v8_dev_results` directory for analysis.

Important: PHOT is development data now, not fresh confirmation. No canonical freeze is created by this run.
