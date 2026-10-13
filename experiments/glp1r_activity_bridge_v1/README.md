# NABU GLP1R activity bridge V1

A bounded experiment testing joint residue/pair estimation against canonical NABU marginal memory equations and same-information supervised baselines. This is measured GLP1R potency prediction, not de novo molecule design. See `PROTOCOL.md` for the frozen hypotheses, exact information contract, and advancement gate.

## Source and evidence

- Puszkarska et al., Nature Chemistry 2024: https://doi.org/10.1038/s41557-024-01532-x
- Official data and code: https://github.com/amp91/PeptideModels, commit `0a23c0bb86ed251acb671db09295aeac7c9c11cf`.
- Unmodified source workbooks: `data/raw/`, with source hashes in `protocol.json` and license in `third_party/`.
- 125 development activity measurements, 19 external records representing 18 primary sequences.
- `DATA_AUDIT.json`, prepared data, and split manifest preserve every source row. The representation reconstructs every sequence exactly, including flanks omitted by the 30-column source alignment.
- `INTEGRITY_TESTS.json` documents pre-outcome synthetic and identity-only checks.
- After execution, `results/REPORT.md` and `results/DECISION.json` contain the outcome; all arms and per-case outputs remain reviewable.

## Reproduce on CPU

Use Python 3.12 and the pinned dependencies. No GPU or paid service is needed. Run from the repository root:

```bash
python -m pip install -r experiments/glp1r_activity_bridge_v1/requirements.txt
python experiments/glp1r_activity_bridge_v1/test_integrity.py
python experiments/glp1r_activity_bridge_v1/predict.py --output experiments/glp1r_activity_bridge_v1/replay
python experiments/glp1r_activity_bridge_v1/predict.py --output experiments/glp1r_activity_bridge_v1/replay_check
python experiments/glp1r_activity_bridge_v1/check_replay.py --original experiments/glp1r_activity_bridge_v1/replay --replay experiments/glp1r_activity_bridge_v1/replay_check
python experiments/glp1r_activity_bridge_v1/evaluate.py --output experiments/glp1r_activity_bridge_v1/replay
```

The committed prepared files are used directly. `prepare.py` is the provenance-preserving extraction recipe; it refuses to overwrite the original prepared state once the original prediction freeze exists. Prediction requires tracked, committed experiment source and data. Both stages refuse to overwrite a run. Choose another fresh `replay_*` directory for an exact repetition.

For the initial execution only, the order is prepare, integrity tests, **protocol/code/data commit**, predict, deterministic replay, compare prediction files, **prediction freeze commit**, evaluate, result commit. Source and data SHA-256 checks guard each boundary. Repeating public retrospective evaluation is reproducibility, not another independent success.

## Interpretation

The candidate is a standard categorical ANOVA kernel ridge model used as a falsifiable restructuring test; no algorithmic novelty is claimed. B2/B3 import unchanged canonical equations, but the 44-slot training-consensus adapter and supervised data contract are new. They are not the complete historical V8.3 pipeline. Every baseline receives the same representation and split.

Only an improvement over the development-selected baseline, the additive ablation, and the frozen uncertainty thresholds on both partitions can advance the candidate to independent replication. Any failed criterion means `DO_NOT_ADVANCE`. A better score than the old memory alone is insufficient. The external study is small, public, related to the same peptide family, and not a prospective wet-lab experiment. No post-result tuning is permitted.
