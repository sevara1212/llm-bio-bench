# How the models fail

Examples of the kinds of wrong answers seen in this benchmark, with question IDs so each one can be looked
up in the results files. Accuracy numbers and full tables are in [`methods.md`](methods.md).

- PBMC3k (8 cell types): `results/scores_pbmc.csv`, run 0, 320 questions per model. Question IDs are from
  `data/questions_pbmc.json`. The PBMC3k prompt says "a cluster from human blood".
- Hao 2021 (30 cell types): `results/scores_hao.csv`, the 600 non-noise questions per model, answers mapped
  with the current judge. The Hao prompt says "a cluster from human PBMCs".

---

## Cell-type annotation, PBMC3k (8 types)

### 1. Cell types that are not in PBMCs

Sometimes a model names a cell type that PBMC preparation removes (neutrophils, erythroid cells) or that is
not in blood at all. The PBMC3k prompt said "human blood" rather than "PBMC", and neutrophils and
reticulocytes are in whole blood, so some of these are partly a prompt problem.

| Model | Question | Truth | Answer | Notes |
|---|---|---|---|---|
| Gemini | `pbmc_c1_noise_4_r0_symbol` | CD4 T cells | neutrophil | 4 random genes plus 6 ribosomal genes |
| Gemini | `pbmc_c7_n_genes_10_r0_ensembl` | Megakaryocytes | Neutrophil | platelet genes not recognised as Ensembl IDs |
| GPT | `pbmc_c3_noise_4_r2_symbol` | CD14 Monocytes | neutrophil | LYZ and S100A9 kept; S100A8/9 are also neutrophil genes |
| Claude | `pbmc_c3_noise_4_r0_ensembl` | CD14 Monocytes | Neutrophil (granulocyte lineage) | same S100 confusion, from Ensembl IDs |
| Claude, Gemini | `pbmc_c6_rank_window_21-30_r0_ensembl` | Dendritic cells | Neutrophil | both models on the same question |
| Claude | `pbmc_c1_noise_4_r1_symbol` | CD4 T cells | Erythroid lineage cells | the random gene SLC40A1 (an iron exporter) pulled it to erythroid |
| Gemini | `pbmc_c1_noise_4_r1_symbol` | CD4 T cells | reticulocyte | same question, same mistake |
| GPT | `pbmc_c7_noise_6_r0_symbol` | Megakaryocytes | endothelial cell | SPARC and GNG11 survived the noise |

Most cases come from noisy or ribosome-heavy gene lists, where one random gene gets too much weight.

### 2. Calling the cluster low quality

GPT answered "unknown / low-quality" 27 times, 25 of them on the CD4 T cluster; Claude and Gemini never
did. Example: `pbmc_c1_ribo_filter_off_r0_symbol` (RPS12, LDHB, RPS25, RPS27, RPS6, ...) -> "Unknown/low-quality
cell". This is a defensible answer, because the CD4 cluster's top markers are ribosomal genes. With the
ribosomal filter on, GPT named the CD4 T cells correctly. Claude and Gemini guessed instead ("T lymphocyte
(naive/resting)", "Naive CD4+ T cell").

### 3. Right lineage, not specific enough

"T cell" for CD4 or CD8 T cells, "monocyte" for CD14 or FCGR3A monocytes. These get half credit under
lenient scoring.

| Model | Strict | Lenient | Share of answers that were right lineage only |
|---|---|---|---|
| Claude | 78.8% | 85.3% | 13.1% |
| Gemini | 72.8% | 77.3% | 9.1% |
| GPT | 71.2% | 75.9% | 9.4% |

Claude gives the most of these answers, so its lead grows under lenient scoring.

### 4. Confusing neighbouring cell types

| Truth | Answer | Count (3 models) |
|---|---|---|
| CD8 T cells | NK cells | 18 |
| FCGR3A Monocytes | Dendritic cells | 14 |
| Dendritic cells | Monocyte (unspecified) | 12 |
| Dendritic cells | B cells | 6 |
| NK cells | CD8 T cells | 5 |

CD8 T and NK cells are confused in both directions. They share cytotoxic genes (GZMA, NKG7, CCL5, CST7),
so this reflects the gene lists more than a model bias. Dendritic cells and B cells share HLA class II genes.

### 5. No usable answer

7 Gemini replies (4 API errors, 3 that hit the 4,000-token limit while reasoning), e.g.
`pbmc_c0_rank_window_31-40_r0_ensembl`. They count as wrong.

### 6. Different answers to the same question

Each PBMC3k question was asked 3 times.

| Model | Same answer in all 3 runs | Right in some runs, wrong in others |
|---|---|---|
| Claude | 83.4% | 11.6% |
| GPT | 81.9% | 10.9% |
| Gemini | 73.8% | 21.2% |

Overall accuracy moved by at most 1.2 points between runs, but about one in five of Gemini's individual
answers changed between correct and wrong. Claude is the most consistent even though it runs without
`temperature=0` (Sonnet 5 does not accept it).

---

## Cell-type annotation, Hao 2021 (30 types)

### 7. Subtypes within a lineage

This is the main failure on the fine types. Strict accuracy by lineage:

| Model | B | CD4 T | CD8 T | DC | Mono | NK | other T (MAIT, gdT, dnT) | other |
|---|---|---|---|---|---|---|---|---|
| Claude | 44% | 2% | 8% | 39% | 78% | 38% | 20% | 70% |
| Gemini | 34% | 8% | 4% | 34% | 65% | 35% | 25% | 64% |
| GPT | 38% | 4% | 6% | 42% | 52% | 33% | 15% | 48% |

No model got any of these 8 types right on any of its 20 questions: B intermediate, B memory, CD4 CTL, CD4
Proliferating, CD4 TEM, CD8 TCM, Treg, dnT.

Most wrong answers name the right lineage but no subtype (B memory -> "B cell (subtype unclear)" 58 times
across the three models; CD8 TEM -> "CD8 T cell (subtype unclear)" 47). The most common wrong specific
answers:

| Truth | Answer | Count (3 models) |
|---|---|---|
| NK_CD56bright | NK | 48 |
| ASDC | pDC | 26 |
| CD4 CTL | NK | 25 |
| CD4 CTL | CD8 TEM | 18 |
| Treg | CD4 Naive | 16 |
| MAIT | CD8 TEM | 13 |

Cytotoxic gene lists (GZMH, GZMK, NKG7, CCL5) tend to become NK or CD8 effector memory, even for CD4 CTLs
and MAIT cells. The models know the cytotoxic programme but not the markers that separate these types
(e.g. KLRB1 and SLC4A10 for MAIT).

### 8. Proliferation recognised, lineage not

For the CD4, CD8 and NK Proliferating types, the models usually answer "proliferating lymphocytes" or
"cycling T/NK cells" without committing to a lineage: Claude did so on 52, Gemini on 50 and GPT on 45 of
the 60 questions. They see the cell-cycle genes (MKI67, STMN1, HMGB2) but not the lineage. The judge maps
these to "Proliferating lymphocyte (lineage unclear)", which earns lenient credit only.

### 9. Cell types that are not in PBMCs (again)

With the prompt saying "PBMCs", this became rare: 6 of Claude's 1,200 Hao answers (all at the highest
noise levels, e.g. neutrophil, basophil), 1 of Gemini's (basophil) and none of GPT's.

---

## Marker lookup baselines (no LLM)

### 10. Housekeeping genes in the database

Some CellMarker lists include ribosomal genes, so a ribosome-only gene list can match them.
`hao_cCD4-Naive_rank_window_1-10_r0_symbol` (RPL30, RPS13, RPL21, RPS27, RPL32, ...) -> Plasmablast,
matched on 6 ribosomal genes. The LLMs at least recognise such lists as uninformative.

### 11. Ties

The top overlap count is tied between two or more labels in 38.7% of CellMarker questions (23.3% for
PanglaoDB), so the tie rule matters. With fractional credit for ties instead of the tie-break, CellMarker's
strict accuracy is 27.5% instead of 28.7%.

---

## Agents

### 12. Reading weak evidence as "no effect" (Task 2)

Given the co-expression of the two genes in control cells as a tool, the agents treated a weak correlation
(|r| <= 0.05) as evidence of no change. This explains 105 of the Claude agent's 122 wrong answers (GPT 81 of
113, Gemini 98 of 121). Weak correlations were just as common for genes that did change as for genes that
did not.

### 13. Merging similar cells when clustering (Task 3)

On PBMC3k the most common error was in clustering: CD8 T cells partly merged into a CD4-majority cluster
(10 of 21 runs), or CD8 T and NK cells merged into one cluster that was then named NK (5 runs). On Hao the
same confusion is larger: 52-62% of CD8 T cells got CD4-lineage labels, across the 15- and 30-step runs.

### 14. Running out of steps (Task 3)

Gemini assumed that Python variables persisted between calls, although the tool description said they do
not, and spent most steps re-running the pipeline after `NameError`s. It used all 15 steps in every run. One
run saved no labels and still ended with a message saying the analysis was complete. On Hao, Claude ran out
of steps in 4 of 7 runs before saving its labels.
