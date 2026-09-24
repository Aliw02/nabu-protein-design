# Historical Experiments and Validation Lineage (NABU V1 - V8.3)

This directory catalog documents the chronological progression and preregistered validation protocols of the project from Version 1 through Version 8.3.

---

## Experiment Lineage Catalog

### 1. `lightning_protein_poc_v1`
* **Objective**: Initial proof-of-concept for structured state memory on single and double mutational scanning data.

### 2. `lightning_protein_poc_ablation_v1`
* **Objective**: Formal ablation study evaluating component-level contributions and anti-leakage audit.

### 3. `lightning_protein_replication_v2`
* **Objective**: Replication study on ion-channel and RNA-binding proteins (KCNH2 and PABP).

### 4. `lightning_gb1_blind_design_v3`
* **Objective**: Prospective blind design on the 4-site combinatorially complete GB1 fitness landscape (149,361 variants).

### 5. `lightning_combinatorial_blind_v3`
* **Objective**: Formal preregistration protocol for cross-assay combinatorial evaluations.

### 6. `lightning_protein_design_v3`
* **Objective**: Large-scale prospective validation across four distinct combinatorial fitness landscapes (GB1, TrpB, PhoQ, TEV).

### 7. `lightning_trpb_blind_design_v4`
* **Objective**: Blind combinatorial candidate selection on the Tryptophan Synthase beta subunit landscape.

### 8. `lightning_protein_design_v4_trpb`
* **Objective**: Data mirror addenda and prospective candidate generation for TrpB.

### 9. `lightning_combinatorial_blind_v4`
* **Objective**: Preregistration and execution framework for the Phototropin (PHOT_CHLRE) 2-5 mutation landscape.

### 10. `lightning_protein_design_v5_phoq`
* **Objective**: Prospective design protocol for the PhoQ sensor kinase signaling domain.

### 11. `lightning_tev_blind_design_v5`
* **Objective**: Prospective design protocol for the Tobacco Etch Virus (TEV) catalytic protease.

### 12. `lightning_clean_blind_confirm_v6`
* **Objective**: Clean blind confirmation protocol on T7 RNA polymerase and Dihydrofolate reductase (DHFR) landscapes.

---

## Modern Epistatic Development Lineage (V8.1 - V8.3)

### 13. `experiments/v8_1_crossfit_higher_order_dev`
* **Objective**: Implementation of cross-fitted triplet and quartet epistatic residual hierarchies with Empirical Bayesian shrinkage.
* **Key Finding**: Demonstrates that higher-order epistasis improves global rank in Phototropin and elite discovery in GB1, while revealing the need to prevent rank distortion on sparse landscapes.

### 14. `experiments/v8_2_objective_tournament`
* **Objective**: 13-arm experimental tournament testing dosage shrinkage, quadratic confidence gates, elite objective splits, rank fusions, and automated out-of-fold landscape gating.
* **Key Finding**: The automated `T_OOF_OBJECTIVE_GATE` achieves undisputed first place across multi-landscape benchmarks (Spearman $\Delta = +0.0065$, 100% positive datasets for both global rank and elite selection).

### 15. `experiments/v8_3_multilandscape_validation`
* **Objective**: Standardized multi-protein validation of the **NABU V8.3 Adaptive Dual-Objective Router** across 4 major benchmark assays:
  - `GB1_V83.csv` (149,361 variants) $\rightarrow$ `RANK_PRESERVING_B3_TOP20_B5_RERANK` (50/50 Top 1% hits, +0.0013 Spearman).
  - `TRPB_V83.csv` (159,129 variants) $\rightarrow$ `GLOBAL_HIGHER_ORDER` (+0.0094 Spearman, Top-1 Regret drops to 0.0341).
  - `PHOQ_V83.csv` (140,517 variants) $\rightarrow$ `B3_PROTECTED_NO_HIGHER_ORDER` (Zero regression safety net).
  - `PHOT_CHLRE_Chen_2023.csv` (14,146 variants) $\rightarrow$ `GLOBAL_HIGHER_ORDER` (+0.0117 Spearman, Top-50 true mean +0.0178).

### 16. `experiments/v8_3_creilov_sealed`
* **Objective**: Preregistered, sealed prospective candidate extrapolation test on the CreiLOV photoreceptor landscape.
* **Protocol**: Trained strictly on low-order mutants ($n \le 3$, 1,175 variants) to predict 13,168 hidden 4-5 mutation variants ($n \in \{4, 5\}$).
* **Key Finding**: Router engaged `B3_PROTECTED` safety net due to low-order training data, achieving $\rho = 0.8968$ (vs $-0.1599$ for shuffled control) with Top-1 Regret of `0.0333` (99.13th percentile).

### 17. `experiments/v8_3_eqfp611_sealed`
* **Objective**: One-shot preregistered higher-order sealed candidate test on the complete 13-site combinatorial fluorescence landscape of eqFP611.
* **Protocol**: Stage A qualification verified 286 triplets and 715 quartets before truth reveal.
* **Outcome**: **`SEALED_CANDIDATE_TEST_PASS`** (Strict improvements on 4 out of 5 primary endpoints: Spearman $\rho$ rose to `0.7894` (+2.93%), Top-10 percentile rose to `91.84%` (+1.16%), Top-1 Regret dropped to `0.1109`).
