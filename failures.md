# Failure types

Running log of *how* models get cell types wrong, with real examples.
Source: `results/scores.csv` (PBMC3k, 320 questions × 3 models, run 0). IDs are question IDs from `data/questions.json`.

Counts are from run 0 only; update after the 3× repeats.

---

## 1. Answer impossible for the sample type

The model names a cell type that PBMC preparation removes (neutrophils, erythroid cells) or that isn't in blood at all.

**Caveat — partly our fault:** the prompt says "a cluster from human blood", **not** "PBMC". Neutrophils and reticulocytes *are* in whole blood, so strictly the model isn't contradicting what it was told. Next prompt version should say "PBMC" and then re-test: if models still say neutrophil, it's a real context-ignoring failure.

| Model | Question | Truth | Predicted | Notes |
|---|---|---|---|---|
| gemini-3.8-flash | `pbmc_c1_noise_4_r0_symbol` | CD4 T cells | neutrophil | 4 random genes + 6 ribosomal genes |
| gemini-3.8-flash | `pbmc_c7_n_genes_10_r0_ensembl` | Megakaryocytes | Neutrophil | top-10 as Ensembl IDs — platelet genes not recognised |
| gpt-5.6-terra | `pbmc_c3_noise_4_r2_symbol` | CD14 Monocytes | neutrophil | LYZ, S100A9 kept — S100A8/9 are shared with neutrophils |
| claude-sonnet-5 | `pbmc_c3_noise_4_r0_ensembl` | CD14 Monocytes | Neutrophil (granulocyte lineage) | same S100 confusion, via Ensembl |
| claude-sonnet-5, gemini | `pbmc_c6_rank_window_21-30_r0_ensembl` | Dendritic cells | Neutrophil | both models, same question |
| claude-sonnet-5 | `pbmc_c1_noise_4_r1_symbol` | CD4 T cells | Erythroid lineage cells (erythroblasts/reticulocytes) | random gene SLC40A1 (iron exporter) pulled it to erythroid |
| gemini-3.8-flash | `pbmc_c1_noise_4_r1_symbol` | CD4 T cells | reticulocyte | same question as above, same mistake |
| gpt-5.6-terra | `pbmc_c7_noise_6_r0_symbol` | Megakaryocytes | endothelial cell | SPARC + GNG11 survive the noise; endothelial cells aren't in blood at all |

**Pattern:** mostly triggered by noise or by a ribosome-heavy list, where one random gene (SLC40A1, SPARC) gets over-weighted.

## 2. Refusing / "low-quality cluster"

The model decides the cluster is technical junk and won't name a cell type.

- **GPT-5.6-terra did this 27 times (25 on the CD4 cluster); Claude and Gemini never did**, e.g. `pbmc_c1_ribo_filter_off_r0_symbol` → "Unknown/low-quality cell" for `RPS12, LDHB, RPS25, RPS27, RPS6, …`.
- Arguably a *reasonable* answer: a ribosome-dominated top-10 list does look like a low-quality cluster. With ribosomal genes filtered out, GPT gets CD4 right every time.
- Claude and Gemini guess instead (Gemini: "naive CD4+ T cell"; Claude: "naive T cell"). Is refusing better or worse than guessing? Worth a sentence in the paper.

## 3. Too vague ("playing it safe")

The answer is the right lineage but not specific enough: "T cell" for CD4/CD8 T cells, "Monocyte" for CD14 / FCGR3A monocytes. Scored 0.5 under lenient scoring, 0 under strict.

| Model | Strict | Lenient | Vague-answer rate |
|---|---|---|---|
| claude-sonnet-5 | 78.8% | 85.3% | **13.1%** |
| gemini-3.8-flash | 72.8% | 77.3% | 9.1% |
| gpt-5.6-terra | 71.2% | 75.9% | 9.4% |

Claude hedges the most — its lead over the others is bigger under lenient scoring.

## 4. Confusing neighbouring cell types

Not bizarre, just the classic hard distinctions. Most common wrong answers:

| Truth | Predicted instead | Count |
|---|---|---|
| CD8 T cells | NK cells | 18 |
| FCGR3A Monocytes | Dendritic cells | 14 |
| Dendritic cells | Monocyte (unspecified) | 12 |
| Dendritic cells | B cells | 6 (HLA class II genes are shared) |
| NK cells | CD8 T cells | 5 |

These go both ways (CD8 ↔ NK), which points at the genes, not a model bias: GZMA, NKG7, CCL5 and CST7 are in both.

## 5. No usable answer

The reply was empty, truncated, or not JSON. 7 cases, all Gemini (4 API errors, 3 hit the 4000-token limit while reasoning), e.g. `pbmc_c0_rank_window_31-40_r0_ensembl`. Counts as wrong; report separately so it isn't mistaken for a biology error.

## 6. Inconsistent between identical calls

Same question, asked 3 times, different answers. From the PBMC3k 3× repeats:

| Model | Same answer all 3 runs | Right in some runs, wrong in others |
|---|---|---|
| claude-sonnet-5 | 83.4% | 11.6% |
| gpt-5.6-terra | 81.9% | 10.9% |
| gemini-3.8-flash | 73.8% | **21.2%** |

Overall accuracy barely moves between runs (±1 point), so the *averages* are stable — but for Gemini, about 1 in 5 individual answers is a coin flip. Note Claude runs without `temperature=0` (Sonnet 5 rejects it) and is still the most consistent.

---

# Hao 2021 CITE-seq (30 fine types) — Claude Sonnet 5, 1200 questions

Prompt now says "PBMCs". Answers mapped to labels by `judge.py` (Claude Sonnet 5; 30/30 on `test_judge.py`,
~5% debatable calls in a 35-answer spot-check).

Headline: strict accuracy falls from **79% (8 types) to 36% (30 types)**. Lenient: 57%.

## 7. Can't split subtypes of one lineage

The main failure on fine types. Strict accuracy by lineage:

| Lineage | CD4 T | other T (MAIT, gdT, dnT) | CD8 T | DC | B | NK | Mono |
|---|---|---|---|---|---|---|---|
| Strict | **2%** | 18% | 24% | 38% | 44% | 45% | 84% |

**Never correct once** (0/40 each): CD4 TCM, CD4 TEM, CD4 CTL, Treg, CD4 Proliferating, B memory,
B intermediate, dnT, ILC. Most common confusions:

| Truth | Answered | Count |
|---|---|---|
| CD4 CTL | CD8 TEM | 30 |
| MAIT | CD8 TEM | 24 |
| gdT | CD8 TEM | 16 |
| dnT | CD8 TCM / CD8 TEM | 24 |
| ILC | CD4 T (subtype unclear) | 15 |

*(Counts above are from the Claude-only run with the earlier judge. With the tightened judge — see
methods.md — plain "cytotoxic CD8" answers now count as "CD8 T cell (subtype unclear)", not CD8 TEM, so
some of these move to that label; the model still answers CD8 either way.)*

**Pattern:** anything cytotoxic (GZMH, GZMK, NKG7, CCL5) becomes "CD8 effector memory", even CD4 CTLs and
MAIT/γδ T cells. The model knows the cytotoxic program but not the markers that tell these apart
(TRDC/TRGC1 for γδ, KLRB1+SLC4A10 for MAIT).

## 8. Proliferation spotted, lineage not

For CD4/CD8/NK Proliferating clusters the model usually says "proliferating lymphocytes (cycling T/NK)".
It sees the cell-cycle genes (MKI67, STMN1, HMGB2) but won't commit to a lineage. The judge gives these
lenient credit via "Proliferating lymphocyte (lineage unclear)".

## 1 (re-check). Impossible answers persist even when the prompt says PBMC

6/1200 answers (0.5%) still named cells that aren't in PBMCs: neutrophil ×4, basophil ×2 — all at noise
level 4 or 6, e.g. `hao_cTreg_noise_4_r0_symbol` → "Neutrophil (granulocyte)". So type 1 is real, but rare
and only under heavy noise.

## Scoring notes for the methods section
- The judge is a Claude model grading Claude answers. It sees only the answer text (no genes, no truth),
  but hand-check ~50 judgements before publishing.
- Judge v1 was too lenient ("conventional DC" → cDC2); v2 is stricter. 102/622 labels changed between versions.

---

## To do
- [x] Re-count everything after the 3× repeats — accuracy per run varies by ≤1.2 points; see type 6
- [ ] Re-check types 1–5 counts using all 3 runs
- [ ] Add "PBMC" to the prompt for the CITE-seq dataset, re-check failure type 1
- [ ] Look at confidence scores: are impossible answers given with high confidence?
