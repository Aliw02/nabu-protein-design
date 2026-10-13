# Commit and execution provenance

Shell Git push did not have credentials. The connected GitHub application published exact Git tree mirrors of the local commits; local and remote commit identifiers differ because commit metadata differs. File trees were checked to be identical before each publication. The local commit identifiers in execution manifests are preserved, not rewritten to the later publication identifiers.

| Boundary | Local commit | Published GitHub commit | Identical tree |
|---|---|---|---|
| Protocol, implementation, raw/prepared data and integrity tests | `5b37722904c7fa431e58d56a22e65e0085528c96` | `9e7ebaf2acb0ee043c17f39849f5ef390ec27faf` | `0c249d5d78c47e4ad0f9c2776ca9424bdb43f9b6` |
| Complete predictions, champion selection and deterministic replay; before external evaluation | `2de0ff8cf7837dc4f199efff0d617140117c72c4` | `0aec592d3562676093fb8fe8ef605ccb0901be2b` | `cfffb45bb99a001139b259df1b0354a5a37a2c03` |

The protocol was committed locally before any activity-data fitting. Its remote mirror was published after fitting and before external evaluation. The complete prediction freeze was committed both locally and on GitHub before invoking the external evaluator. This is retrospective public-data evaluation with an auditable process boundary, not a globally blind or prospectively registered clinical experiment.

`results/PRE_RUN_MANIFEST.json`, `results/PREDICTION_FREEZE.json`, `results/DETERMINISM.json`, and `results/DECISION.json` record execution times, file hashes, source hashes, deterministic checks, and the frozen decision. The full prediction stage took 21.1403 seconds on the recorded single-thread CPU environment; replay is recorded separately. Model times are descriptive and specific to this tiny dataset and machine.

No model code, protocol, data, grid, metric, or gate was changed after the pre-fit commit. `INTERPRETATION.md` is explicitly post-outcome commentary. Exact reproduction on the published tree checks file hashes, not equality between local and remote commit metadata.
