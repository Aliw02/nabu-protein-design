# NABU Phase-1 Benchmark Record

This document summarizes the final Phase-1 evidence for NABU V8.3.

Phase 1 is closed. Historical raw outputs and preregistrations remain authoritative where a summary and a raw artifact differ.

## 1. Core multi-landscape router evidence

| Landscape | Router behavior | B3 Spearman | V8.3 Spearman | Main observation |
| --- | --- | ---: | ---: | --- |
| PHOT | `GLOBAL_HIGHER_ORDER` | 0.8995 | **0.9112** | Higher-order useful globally |
| TrpB | `GLOBAL_HIGHER_ORDER` | 0.2902 | **0.2996** | Global rank and elite discovery improved |
| GB1 | `RANK_PRESERVING_B3_TOP20_B5_RERANK` | 0.4018 | **0.4030** | Global HO harmful; elite-only rerank useful |
| PhoQ | `B3_PROTECTED_NO_HIGHER_ORDER` | **0.5420** | **0.5420** | Router prevented higher-order regression |

These landscapes established why V8.3 requires routing rather than universal use of B5.

## 2. CreiLOV sealed pairwise extrapolation

Protocol:

- visible: mutation count <= 3,
- hidden: mutation count 4 or 5,
- hidden rows: 13,168,
- higher-order memory at Stage A: no supported triplet/quartet corrections.

Result:

- Router mode: `B3_PROTECTED_NO_HIGHER_ORDER`
- Spearman: **0.8968**
- Top-10 mean percentile: **0.9788**
- Top-50 mean percentile: **0.9513**
- Top-50 true Top-1% hits: **13**
- Normalized Top-1 regret: **0.0333**
- Shuffled control Spearman: **-0.1599**

Interpretation: strong sealed evidence for B3 pairwise extrapolation into unseen 4/5-mutation variants, but not a higher-order validation because higher-order memories were empty.

## 3. eqFP611 sealed higher-order validation

Stage-A qualification:

- visible rows: 1,599
- hidden rows: 689
- triplet entries: **286**
- quartet entries: **715**
- visible OOF B3: 0.7266
- visible OOF B4: 0.7613
- router: `GLOBAL_HIGHER_ORDER`

Hidden evaluation:

| Endpoint | B3 | V8.3 | Result |
| --- | ---: | ---: | --- |
| Spearman | 0.7601 | **0.7894** | strict gain |
| Top-10 mean percentile | 0.9068 | **0.9184** | strict gain |
| Top-50 mean percentile | 0.8983 | **0.9015** | strict gain |
| Top-50 true Top-1% hits | 3 | 3 | no regression |
| Normalized Top-1 regret | 0.1281 | **0.1109** | strict gain |

Final status: `SEALED_CANDIDATE_TEST_PASS`.

Shuffled-label Spearman: **0.0514**.

## 4. RhlA sample-efficiency benchmark

Fixed hidden triple-mutant test: **540 variants**.

Passive V8.3:

| Training budget | Rows | Spearman | Top-10 mean pct | Top-50 mean pct |
| --- | ---: | ---: | ---: | ---: |
| 5% | 118 | 0.5685 | 0.8391 | 0.7019 |
| 10% | 235 | 0.6785 | 0.7787 | 0.7660 |
| 20% | 469 | 0.6900 | 0.8400 | 0.7607 |
| 40% | 938 | **0.7620** | **0.8422** | **0.8125** |

The preregistered practical threshold was first satisfied by **40%**, not 20%.

Negative controls:

- random-hash Spearman: approximately -0.024
- shuffled-label Spearman: approximately 0.119

## 5. RhlA Active Acquisition development

This experiment was **retrospective controller development** because RhlA truth had already been revealed before the controller was designed.

The 50/50 exploration/exploitation controller produced:

- Active 20% Spearman: **0.7202**
- Active 20% Top-10 mean percentile: **0.8907**
- Active 20% Top-50 mean percentile: **0.8254**
- Active 20% Top-50 Top-1% hits: **3**

Active 20% reached the frozen Passive 40% practical threshold:

- measurement reduction: **50%**
- measurement efficiency factor: **2.0×**

This was useful development evidence, not fresh prospective proof.

## 6. CR9114-H1 final fresh Phase-1 gate

This was the final one-shot fresh Phase-1 experiment.

Candidate scope:

- training labels: measured 1–5 mutation variants,
- hidden evaluation: 4/5-mutation variants,
- fixed hidden rows: **1,673**,
- common scoreability coverage: **100%**.

### Higher-order qualification before reveal

- Passive 40 triplets: **546**
- Passive 40 quartets: **1,698**
- Active 20 triplets: **545**
- Active 20 quartets: **1,024**

Qualification passed.

### Core V8.3 final gate

Passive 40%:

| Endpoint | B3 | V8.3 |
| --- | ---: | ---: |
| Spearman | 0.9271 | **0.9387** |
| Top-10 mean percentile | 0.9179 | **0.9297** |
| Top-50 mean percentile | 0.8902 | **0.9047** |
| Top-50 true Top-1% hits | 4 | **5** |
| Normalized Top-1 regret | 0.0476 | 0.0476 |

Core gate:

- all-no-regression: **true**
- strict gains: **4/5**
- core pass: **true**

Negative controls:

- random-hash Spearman: **-0.0285**
- shuffled-label V8.3 Spearman: **-0.1111**

### Active controller final gate

Against Passive 40% V8.3:

- Active 5%: failed threshold
- **Active 10%: passed threshold**
- Active 20%: failed threshold
- Active 40%: failed threshold

Active 10% retained **91.19%** of Passive 40% Spearman while satisfying all frozen Top-10, Top-50, and Top-1%-hit gap criteria.

This corresponds to:

- **75% fewer measurements**
- **4.0× measurement efficiency**

However, the preregistered final target required **Active 20%**, and Active 20% failed. Therefore the formal status is:

`PHASE1_FINAL_GATE_CLOSED_WITH_LIMITATION`

The limitation is specifically the non-monotonic Active Acquisition policy, not the V8.3 core.

## 7. Final Phase-1 interpretation

Established:

- B3 is a strong protected pairwise backbone.
- Cross-fitted higher-order memory can improve unseen combinatorial ranking.
- Visible OOF routing successfully distinguishes global, elite-only, and protected regimes.
- The frozen core passed a fresh high-order 4/5-mutation final gate.
- Active selection can reduce measurements substantially at some budgets.

Not established:

- universal monotonic superiority of the current 50/50 Active Acquisition controller,
- de novo full-protein sequence generation,
- structure-conditioned design,
- wet-lab validation of NABU-proposed novel variants outside published benchmark landscapes.

Phase 2 starts from this boundary.
