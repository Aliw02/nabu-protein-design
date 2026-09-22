# Historical Experiments and Validation Lineage

This directory catalog documents the chronological progression and preregistered validation protocols of the project from Version 1 through Version 6.

---

## Experiment Directory Catalog

### 1. `lightning_protein_poc_v1`
* **Objective**: Initial proof-of-concept for structured state memory on single and double mutational scanning data.
* **Key Artifacts**:
  * `summary.json`: Initial validation metrics.
  * `predictions_revealed.csv`: Unblinded prediction scores.
  * `metrics.csv`: Primary rank correlation summary.

### 2. `lightning_protein_poc_ablation_v1`
* **Objective**: Formal ablation study evaluating component-level contributions and anti-leakage audit.
* **Key Artifacts**:
  * `run_ablation.py`: Full ablation execution pipeline.
  * `ablation_metrics.csv`: Quantitative impact of individual algorithmic layers.
  * `component_deltas.csv`: Residual delta contributions.
  * `leakage_and_reproduction_audit.json`: Deterministic data leakage verification.

### 3. `lightning_protein_replication_v2`
* **Objective**: Replication study on ion-channel and RNA-binding proteins (KCNH2 and PABP).
* **Key Artifacts**:
  * `preregistration.json`: Protocol freeze prior to evaluation.
  * `replication_summary.json`: Consolidated multi-protein results.
  * `kcnh2_corrected_group_metrics.json`: High-throughput assay validation.

### 4. `lightning_gb1_blind_design_v3`
* **Objective**: Prospective blind design on the 4-site combinatorially complete GB1 fitness landscape (149,361 variants).
* **Key Artifacts**:
  * `PREREGISTRATION_PRE_REVEAL.json`: Cryptographic protocol commitment.
  * `CANDIDATE_FREEZE_PRE_REVEAL.json`: SHA-256 frozen candidate manifest.
  * `SHUFFLE_CONTROL_FREEZE_PRE_REVEAL.json`: Randomized control freeze.
  * `results_summary.json`: Unblinded hit rate analysis.

### 5. `lightning_combinatorial_blind_v3`
* **Objective**: Formal preregistration protocol for cross-assay combinatorial evaluations.
* **Key Artifacts**:
  * `preregistration.json`: Standardized split rules (FNV1a hash mod 10) and pass gates.

### 6. `lightning_protein_design_v3`
* **Objective**: Large-scale prospective validation across four distinct combinatorial fitness landscapes (GB1, TrpB, PhoQ, TEV).
* **Key Artifacts**:
  * `PREREGISTRATION_V3.json`: Primary multi-protein preregistration document.
  * `V3_RESULTS.json`: Consolidated quantitative metrics across all 547,161 variants.
  * `V3_METRICS.csv`: Per-landscape candidate breakdown.
  * `V3_SHAREABLE_REPORT.md`: Comprehensive evaluation report.

### 7. `lightning_trpb_blind_design_v4`
* **Objective**: Blind combinatorial candidate selection on the Tryptophan Synthase beta subunit landscape.
* **Key Artifacts**:
  * `PREREGISTRATION_PRE_REVEAL.json`: TrpB-specific candidate generation protocol.
  * `CANDIDATE_FREEZE_PRE_REVEAL.json`: Frozen candidate manifest with cryptographic hashes.
  * `visible_training.csv`: Non-overlapping visible training windows.

### 8. `lightning_protein_design_v4_trpb`
* **Objective**: Data mirror addenda and prospective candidate generation for TrpB.
* **Key Artifacts**:
  * `PREREGISTRATION_V4_TRPB.json`: Extended protocol definition.
  * `DATA_MIRROR_ADDENDUM.json`: Source tracking and integrity audit.

### 9. `lightning_combinatorial_blind_v4`
* **Objective**: Preregistration and execution framework for the Phototropin (PHOT_CHLRE) 2-5 mutation landscape.
* **Key Artifacts**:
  * `preregistration.json`: Protocol definition and evaluation thresholds.

### 10. `lightning_protein_design_v5_phoq`
* **Objective**: Prospective design protocol for the PhoQ sensor kinase signaling domain.
* **Key Artifacts**:
  * `PREREGISTRATION_V5_PHOQ.json`: Assay parameters and scoring definitions.
  * `VISIBLE_COVERAGE_ADDENDUM.json`: Positional coverage constraints.

### 11. `lightning_tev_blind_design_v5`
* **Objective**: Prospective design protocol for the Tobacco Etch Virus (TEV) catalytic protease.
* **Key Artifacts**:
  * `PREREGISTRATION_PRE_REVEAL.json`: Catalytic activity thresholds and candidate budgets.

### 12. `lightning_clean_blind_confirm_v6`
* **Objective**: Clean blind confirmation protocol on T7 RNA polymerase and Dihydrofolate reductase (DHFR) landscapes.
* **Key Artifacts**:
  * `PREREGISTRATION_PRE_SEARCH.json`: Pre-search cryptographic commitment for T7 and DHFR benchmarks.
