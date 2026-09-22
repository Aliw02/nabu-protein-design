# Lightning Architecture — Blind Combinatorial Protein Design POC V3

## Result

A frozen Lightning-derived protein adapter generated and ranked **unseen 4-site amino-acid combinations** from visible experimental measurements on four independent combinatorial fitness landscapes: TrpB, GB1, PhoQ, and TEV.

Across the first 10 proposed candidates per landscape (40 total):

| Method | Hidden top-1% hits | Hit rate | Mean hidden percentile | Worst-landscape mean percentile |
|---|---:|---:|---:|---:|
| B3 raw pair interaction | 34/40 | 85.0% | 99.31% | 98.49% |
| B5 Lightning consensus | 33/40 | 82.5% | 99.46% | 99.18% |
| B2 additive memory | 23/40 | 57.5% | 98.40% | 97.61% |
| Hamming-1 neighbor mean | 15/40 | 37.5% | 79.32% | 44.10% |
| Random hash | 1/40 | 2.5% | 56.32% | 44.68% |
| Shuffled-fitness B5 | 0/40 | 0.0% | 51.85% | 38.08% |

The full B5 arm therefore placed **33/40 = 82.5%** of its first ten hidden designs in the top 1% of the experimentally measured hidden landscape. The raw pair-interaction arm placed **34/40 = 85.0%**. The shuffled-label control placed **0/40**, and the deterministic random control placed **1/40**.

For B5, the preregistered 1,000-sample random-null tests gave one-sided p≈0.001 for the Top-10 mean and Top-1% hit count on **all four landscapes**. Best-candidate random-null p-values were ≈0.002 (TrpB), ≈0.002 (GB1), ≈0.010 (PhoQ), and ≈0.029 (TEV).

## What was actually tested

The predictor was not given the hidden candidate fitness values. Visible measurements were converted into:

- exact four-site state identities;
- visible-only site-state memory;
- visible-only pairwise residual interactions;
- relative within-site and within-pair rankings;
- a B5 consensus score that rewards agreement across four memory and six resonance views.

Every arm generated a ranked list from the full 20^4 candidate space, excluding visible variants. Ranked candidate lists were frozen with SHA-256 hashes before hidden outcomes were evaluated. During reveal, every candidate-list hash was recomputed and required to match the frozen hash.

## Independent landscapes

- **TrpB**: 16,000 visible / 143,129 hidden measured variants.
- **GB1**: 15,000 visible / 134,361 hidden measured variants.
- **PhoQ**: 14,068 visible / 126,449 hidden measured variants.
- **TEV**: 15,910 visible / 143,222 hidden measured variants.

Total hidden measured variants evaluated: **547,161**.

## Hard controls

**Shuffled fitness:** Visible genotype-to-fitness correspondence was deterministically destroyed and the whole B5 pipeline was rerun. It produced **zero top-1% hits in the first 10 candidates across all four landscapes**.

**Random ranking:** A fixed FNV1a genotype hash produced only **1/40** top-1% hits.

**Hamming-1 local baseline:** A strong local-neighbor baseline reached **15/40** top-1% hits, substantially below B3/B5 at the strict Top-10 budget.

## The important negative result

The experiment does **not** support the claim that the current full B5 consensus is uniquely better than a simpler pairwise model.

B3 raw pair interaction beat B5 on 3/4 landscapes by Top-10 mean percentile and was generally stronger at K=50 and K=100. B5 had a slightly higher cross-landscape average Top-10 percentile and a better worst-landscape Top-10 percentile, but this does not erase B3's broader advantage.

That is the main architectural diagnosis from V3: **the protein-transfer signal is real, but most of the transferable design power currently sits in structured state memory + pairwise interaction/resonance. The present B5 consensus wrapper is not yet justified as the superior final selector.**

## Shareable claim supported by these results

> A frozen structured-memory and interaction-based reasoning pipeline prospectively generated unseen combinatorial protein variants that were strongly enriched for experimentally high-fitness sequences across four independent protein landscapes. In the first ten proposals per landscape, the full consensus arm placed 82.5% of candidates in the hidden top 1%, versus 37.5% for a Hamming-1 local baseline, 2.5% for deterministic random selection, and 0% after shuffling visible fitness labels. A simpler pairwise ablation performed slightly better overall (85.0%), identifying the current pairwise interaction layer as the principal transferable mechanism.

## Claim boundary

This is **computational retrospective/prospective-style validation against already-existing experimental measurements**, not new wet-lab validation. It demonstrates unseen combinatorial candidate generation and blind outcome recovery under frozen candidate hashes. It does not establish that arbitrary de-novo proteins outside these measured landscapes will function experimentally.

## Integrity artifacts

- Main V3 preregistration: `V3_PREREGISTRATION.json`
- PhoQ preregistration and source amendment: `V3_PHOQ_PREREGISTRATION.json`, `V3_PHOQ_SOURCE_MIRROR_AMENDMENT.json`
- TEV preregistration: `V3_TEV_PREREGISTRATION.json`
- Pre-reveal freezes: `TRPB_V3_PRE_REVEAL_FREEZE.json`, `GB1_V3_PRE_REVEAL_FREEZE.json`, `PHOQ_V3_PRE_REVEAL_FREEZE.json`, `TEV_V3_PRE_REVEAL_FREEZE.json`
- Complete metric bundle: `V3_RESULTS.json`, `V3_METRICS.csv`
