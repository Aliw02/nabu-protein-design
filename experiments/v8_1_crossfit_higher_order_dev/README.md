# NABU V8.1 — Cross-Fitted Higher-Order Development

This branch tests whether V8's higher-order gain survives removal of visible-row self-influence from residual targets.

## PHOT development run

From the repository root:

```powershell
py -m pip install pandas numpy scipy
py experiments/v8_1_crossfit_higher_order_dev/run_v8_1_crossfit_dev.py PHOT_CHLRE_Chen_2023.csv --out nabu_v8_1_crossfit_dev_results
```

Expected outputs:

```text
nabu_v8_1_crossfit_dev_results/
  V8_1_CROSSFIT_DEV_RESULTS.json
  V8_1_CROSSFIT_DIAGNOSTICS.json
  V8_1_CROSSFIT_REVEALED_PREDICTIONS.csv
  V8_1_CROSSFIT_TOP50.csv
```

This is development evidence only. Do not freeze or retune from the PHOT result alone.

The same runner accepts another benign combinatorial CSV without code changes.
