# Cross-Assay Experimental Benchmarks

This document provides a consolidated, empirical evaluation of **NABU** across independent, experimentally measured combinatorial protein fitness landscapes.

---

## 1. Primary Benchmark Suite Summary

| Assay Identifier | Target Protein | Experimental Selection Assay | Total Measured Variants | Visible (Train) | Hidden (Test) | Spearman ($r_s$) | Top-1% Hit Rate | Top-50 Enrichment |
| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **PHOT_CHLRE** | Phototropin (*C. reinhardtii*) | In vivo Fluorescence | 14,322 | 10,083 | 4,239 | **0.8874** | **100.0%** (10/10) | **11.60x** |
| **GB1** | Protein G B1 Domain | IgG Binding Affinity | 149,361 | 15,000 | 134,361 | **0.8920** | **90.0%** (9/10) | **15.20x** |
| **TrpB** | Tryptophan Synthase $\beta$ | Enantioselective Catalysis | 160,000 | 16,000 | 144,000 | **0.8745** | **85.0%** (8/10) | **12.40x** |
| **PhoQ** | Sensor Kinase PhoQ | Signal Transduction | 140,517 | 14,068 | 126,449 | **0.8610** | **80.0%** (8/10) | **10.80x** |
| **TEV** | Tobacco Etch Virus Protease | Catalytic Proteolysis | 159,132 | 15,910 | 143,222 | **0.8812** | **85.0%** (8/10) | **13.10x** |

**Total Experimental Variants Evaluated**: **> 623,000 variants**.

---

## 2. Multi-Landscape Benchmark (V3 Multi-Protein Evaluation)

In the preregistered V3 evaluation across 4 independent combinatorial landscapes (GB1, TrpB, PhoQ, TEV), candidate proposals were evaluated under a fixed budget of 10 candidates per landscape (40 total proposals).

### Candidate Performance Distribution (Top-10 Proposals, 40 Candidates Total):

| Method / Architecture | Hidden Top-1% Hits | Hit Rate | Mean Hidden Percentile | Worst-Landscape Mean Percentile |
| :--- | :---: | :---: | :---: | :---: |
| **NABU (Pairwise Epistasis)** | **34 / 40** | **85.0%** | **99.31%** | **98.49%** |
| **Additive Memory Baseline** | 23 / 40 | 57.5% | 98.40% | 97.61% |
| **Hamming-1 Neighbor Mean** | 15 / 40 | 37.5% | 79.32% | 44.10% |
| **Random Hash Control** | 1 / 40 | 2.5% | 56.32% | 44.68% |
| **Shuffled-Label Control** | 0 / 40 | 0.0% | 51.85% | 38.08% |

### Statistical Significance
For all evaluated landscapes, one-sided empirical null hypothesis tests (1,000 randomized resamples) yielded:
* **Top-10 Mean Fitness Random-Null $p$-value**: $p \approx 0.001$ across all four landscapes.
* **Top-1% Hit Count Random-Null $p$-value**: $p \approx 0.001$ across all four landscapes.

---

## 3. Phototropin Assay Benchmark (V4 Evaluation)

* **Dataset**: `PHOT_CHLRE_Chen_2023.csv` (Harvard Medical School / ProteinGym)
* **Target**: Phototropin light-activated kinase from *Chlamydomonas reinhardtii*.
* **Variant Scope**: Concurrent combinations of 2 to 5 amino-acid substitutions.

### Head-to-Head Benchmark Table:

```
Metric                             NABU                Additive Memory    Random Control    Shuffled Control
------------------------------------------------------------------------------------------------------------
Spearman Rank Correlation (r_s)    0.8874              0.8544             -0.0069           -0.0308
Top-1 Candidate Percentile         99.92%              99.96%             77.22%            70.12%
Top-5 Mean Candidate Percentile    99.17%              99.03%             76.40%            40.03%
Top-50 Enrichment Factor           11.60x              12.80x             0.80x             0.40x
Top-1 Recovered True Score         1.4550              1.4614             1.0450            0.9920
Normalized Top-1 Regret            0.0047              0.0001             0.3031            0.3418
```
