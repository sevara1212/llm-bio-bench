# Methods notes

Working notes for the methods section. Numbers are from `results/scores_hao.csv` unless stated.

## Task 1 pilot: PBMC3k (8 cell types)

- **Data** (`src/prepare_data.py`): `scanpy.datasets.pbmc3k()`, the scanpy tutorial pipeline (cells with
  200-2,500 genes and < 5% mitochondrial counts; highly variable genes by the tutorial's dispersion
  thresholds; counts and mitochondrial share regressed out; PCA; 10 neighbours on 40 PCs; Leiden at
  resolution 1.0), giving the tutorial's 8 cell types (CD4 T, CD8 T,
  B, CD14 monocytes, FCGR3A monocytes, NK, dendritic cells, megakaryocytes). Clusters were labelled from
  canonical tutorial markers (z-scored across clusters) and checked by eye. Top-50 markers per cluster
  (Wilcoxon).
- **Questions** (`src/make_questions.py pbmc`): 8 clusters x 20 question variants (the same four
  difficulty knobs as Hao, incl. noise) x 2 gene formats = 320. Prompt: "...a cluster from human blood"
  (the Hao prompt says "PBMCs"). **Each question was asked 3 times** (runs 0-2): n = 960 answers per
  model.
- **Scoring:** keyword rules (`src/scoring.py`, `RULES`), not the LLM judge - with 8 broad types each
  answer is mapped by the first matching pattern. Strict = exact type; lenient = 0.5 for "T cell" on a
  CD4/CD8 cluster or "monocyte" on a CD14/FCGR3A cluster. The rules read the whole answer text, so they
  do not apply the primary-answer / hedge rules of the Hao judge.
- Models: Claude Sonnet 5 via the Anthropic API; GPT-5.6 Terra and Gemini 3.8 Flash via OpenRouter.
  Random guess: 12.5% (1/8).

| Method (n = 960) | Strict | Lenient | Strict symbol | Strict Ensembl |
|---|---|---|---|---|
| Claude Sonnet 5 | 78.6% | 85.1% | 82.3% | 75.0% |
| Gemini 3.8 Flash | 72.6% | 77.4% | 82.1% | 63.1% |
| GPT-5.6 Terra | 70.5% | 75.2% | 75.4% | 65.6% |
| CellMarker lookup (tie-break) | 64.4% | 70.3% | 64.4% | 64.4% |
| PanglaoDB lookup (tie-break) | 36.9% | 40.9% | 36.9% | 36.9% |
| Random guess | 12.5% | - | - | - |

- Lookups with fractional tie credit: CellMarker 61.0%, PanglaoDB 41.9% strict (they are deterministic,
  so the 3 runs are identical).
- **Run-to-run stability** (strict accuracy moved by at most 1.2 points between runs): per run (0 / 1 / 2): Claude 78.8 / 78.4 / 78.8%, Gemini
  72.8 / 72.2 / 72.8%, GPT 71.2 / 70.0 / 70.3%. Same mapped answer in all 3 runs: Claude 83.4%, GPT
  81.9%, Gemini 73.8% of questions; right in some runs and wrong in others: Claude 11.6%, GPT 10.9%,
  Gemini 21.2%.
- **Difficulty knobs**, strict, both formats and 3 runs pooled (n = 48 per cell; noise levels 2-6 pool
  3 random draws, n = 144):

  | Knob | Level | Claude | Gemini | GPT | CellMarker lookup |
  |---|---|---|---|---|---|
  | marker ranks shown | 1-10 | 88% | 85% | 85% | 62% |
  | | 11-20 | 58% | 48% | 52% | 62% |
  | | 21-30 | 65% | 52% | 52% | 50% |
  | | 31-40 | 65% | 38% | 38% | 12% |
  | genes shown | 1 | 50% | 52% | 48% | 100% |
  | | 3 | 77% | 71% | 62% | 88% |
  | | 5 | 77% | 75% | 69% | 75% |
  | | 10 | 88% | 85% | 79% | 62% |
  | random genes among 10 | 0 | 85% | 85% | 79% | 62% |
  | | 2 | 87% | 88% | 81% | 62% |
  | | 4 | 82% | 70% | 77% | 62% |
  | | 6 | 77% | 70% | 65% | 71% |
  | ribosomal/MT genes removed | off | 85% | 88% | 83% | 62% |
  | | on | 98% | 90% | 94% | 62% |

- **Ribosomal filter and the CD4 cluster:** the CD4 T cluster's top-10 markers are ribosomal genes.
  Lenient score on that cluster without / with the ribosomal/MT filter (n = 6 each: 2 formats x 3
  runs; small): Claude 0.50 / 0.92, Gemini
  1.00 / 0.83, GPT 0.00 / 1.00, CellMarker lookup 1.00 / 1.00.
- Cost (logged): Claude $2.24, Gemini $3.46, GPT $5.04 (960 answers each).

## Task 1 data (Hao 2021)

- CELLxGENE dataset "nygc multimodal pbmc" (Hao et al. 2021, 161,764 cells, `celltype.l2` labels).
  Doublets are dropped, leaving 30 cell types; up to 500 cells per type are subsampled (seed 0) and the
  authors' labels are used directly (no re-clustering). Top-50 markers per type from a Wilcoxon test on
  log-normalised counts (`src/prepare_hao.py`).
- Questions (`src/make_questions.py hao`): 30 types x 20 question variants x 2 gene formats = 1,200;
  the 600 non-noise questions are the main comparison.

## Scoring free-text answers (Hao 2021, 30 cell types)

Each model answer (free text) is mapped to one of the 30 Hao `celltype.l2` labels, a lineage-only
label (e.g. "CD4 T cell (subtype unclear)", "Proliferating lymphocyte (lineage unclear)"),
"Not found in PBMCs", or "Unknown". Mapping is identical for every model:

1. **Exact match** (`judge.rule_match`): after lower-casing, trimming and dropping a plural "s", an answer
   that is exactly a label name or one of a short list of unambiguous synonyms is mapped directly.
   No LLM is involved. (60 of 882 distinct answers.)
2. **LLM judge** for everything else: Claude Sonnet 5 (`claude-sonnet-5`, effort `low`, structured
   output constrained to the label list). The judge sees only the answer text - not the model that
   wrote it, the true label or the genes - so it maps wording and does not grade. Each distinct
   answer string is judged once and the decision is reused for every occurrence, so all models get
   the same label for the same wording.

Strict accuracy = exact fine type. Lenient = 1 for exact, 0.5 for the right lineage (`celltype.l1`)
or a lineage-only answer, 0 otherwise.

Models are compared only on (question, run) pairs answered by all three models: n = 600 per model
(30 cell types x 20 non-noise question variants), n = 30 per point in the difficulty plots.

## Judge validation

- **Test set** (`src/test_judge.py`): hand-written answer strings with known labels. With 34 cases the
  judge scored 34/34, 33/34 and 31/34 on three repeated runs (96% overall); residual misses were almost
  all "proliferating/cycling NK cells" mapped to NK or NK_CD56bright. Extended to 41 cases on 2026-09-28
  for the "likely" rule (below); latest run 40/41, the miss being the open "memory/effector" case.
- **Hand-check, before the fix:** 30 randomly sampled judge decisions (distinct answers, model names
  hidden; `results/hao_judge_check.csv`) were checked by hand. **The judge agreed on 27/30.** The 3
  disagreements:
  1. "CD8+ memory/effector T cells (KLRB1+ GZMK+ subset, possibly MAIT-like)" -> judge said MAIT;
     should be CD8 TEM (the judge credited the hedge, not the primary answer).
  2. "cytotoxic CD8+ T cell" -> judge said CD8 TEM; should be CD8 T cell (subtype unclear).
  3. "CD8+ cytotoxic T cell" -> same as 2.
- **What was changed after the hand-check:**
  1. The judge maps the **primary answer only** and ignores hedges ("possibly X", "maybe X", "or Y",
     "X-like", "e.g. X"). Words describing cell state (proliferating, naive, memory, effector) count as
     part of the primary answer.
  2. **"Cytotoxic CD8" alone = CD8 T cell (subtype unclear).** CD8 TEM only when the answer says
     effector memory, effector or TEM.
  3. The structured output now asks for the primary answer (free text) before the label; if that
     primary answer is itself exactly a label or synonym, it is used. This fixes a label-selection slip
     seen in testing (primary answer "NK Proliferating", label "NK_CD56bright").
  Only cached decisions the new rules could affect were re-judged (561 of 822 distinct answers
  containing "cytotoxic", "effector" or a hedge word); the other 261 were reused. After the fix all 3
  hand-check disagreements are labelled as the hand-check says.
- **"likely" rule (decided 2026-09-28):** "likely X" (or "probably X") counts as the answer, so X is
  judged, e.g. "CD4 T cell (likely naive)" -> CD4 Naive. "possibly X" and "X or Y" / "X/Y" alternatives
  are ignored and only the primary answer is judged; if "likely" introduces alternatives ("likely CD8+
  central/effector memory"), only what they share counts (-> CD8 T cell (subtype unclear)).
  - Check against how the judge had actually scored: it was inconsistent - e.g. "B cells (likely naive
    B cells)" -> B naive (likely counted) but "CD4 T cell (likely naive)" -> CD4 T cell (subtype unclear)
    (likely ignored). All 274 cached answers containing "likely"/"probably" were re-judged under the
    written rule; 33 of 1,114 distinct labels changed.
  - Effect (600 Hao questions): plain Claude strict 32.5% -> 32.7%, lenient 54.9% -> 54.7%; GPT-5.6
    strict 27.2% -> 27.0%, lenient 47.8% -> 47.5%; 96-question subset, plain Claude lenient 53.6% ->
    52.6%. Gemini, the agents and the lookups are unchanged. All numbers in this file use the new rule;
    the "Effect of the fix" table above is historical (the state after the earlier fix).
  - **Not yet done:** answers with "/" or "or" alternatives but no "likely" (509 distinct answers) still
    carry labels from the earlier rule text, which named "or Y" but not "X/Y"; some are inconsistent
    (e.g. three "effector/central memory CD8" variants -> CD8 TCM, one -> CD8 T cell (subtype unclear)).
    Re-judging them was interrupted when the Anthropic credit ran out; their previous labels were
    restored, so every score here is from one consistent cache state.
- **Effect of the fix on accuracy** (same 600 questions per model):

  | Model | Strict before | Strict after | Change | Lenient change | Labels changed | Right/wrong flips |
  |---|---|---|---|---|---|---|
  | Claude Sonnet 5 | 35.3% | 32.5% | -2.8 pts | -1.8 pts | 47 | 17 |
  | Gemini 3.8 Flash | 30.2% | 29.8% | -0.3 pts | -0.2 pts | 2 | 2 |
  | GPT-5.6 Terra | 29.3% | 27.2% | -2.2 pts | -0.9 pts | 42 | 13 |

  Almost all changes are "CD8 TEM" -> "CD8 T cell (subtype unclear)" (51 labels). Claude and GPT often
  answer "cytotoxic CD8 T cell"; Gemini rarely does, so it is barely affected. The model ranking is
  unchanged.
- **Judge determinism:** before the fix, re-judging all of Claude's answers from scratch with unchanged
  instructions changed 3.8% of labels (23/600) but only 0.7% of correct/incorrect outcomes (4/600);
  strict accuracy moved by 0.3 points.

## Baselines (no LLM)

- **Random guess:** uniform over the 30 fine types -> strict 3.3%, lenient 7.8% (PBMC3k: 12.5% strict).
- **Marker lookup** (`src/baseline_marker_lookup.py`), run with two databases:
  - **PanglaoDB** (`PanglaoDB_markers_27_Mar_2020.tsv.gz`, panglaodb.se, downloaded 2026-09-24), human
    entries of organ "Immune system"/"Blood" plus "Hematopoietic stem cells".
  - **CellMarker** (`human_cell_marker.txt`, downloaded by hand 2026-09-24; newer than CellMarker 2.0,
    includes 2025 papers), human, tissue class Blood, disease Normal. **Rows from Hao et al. 2021
    (PMID 34062119) were excluded** - CellMarker contains 2,014 marker rows derived from that very
    dataset (B, CD4 T, CD8 T, DC, monocyte, NK), which would leak the answer key.
  - Method, identical for both: candidate types are database types mapped by hand to the Hao labels
    (PanglaoDB table in the script; CellMarker in `configs/cellmarker_to_hao.csv`, 161 names); "known
    markers" = each type's top 50 genes by evidence (CellMarker: number of distinct papers; PanglaoDB:
    human sensitivity), so large generic lists (CellMarker "B cell": 1,298 genes) can't win by size;
    Ensembl IDs are converted to symbols with the dataset's own ID table; pick the type with the most
    given genes among its markers; no overlap -> Unknown.
  - **Ties** on the top count are broken by the lowest sum of the matched markers' ranks (stronger
    markers win), then alphabetically. Because the last step is arbitrary, the lookups are also scored
    with **fractional credit**: each of the k distinct labels tied on the top count gets 1/k (several
    database names can map to one label, so k counts labels, not names). Both versions are reported.
  - Types a database can't split (e.g. PanglaoDB "Monocytes") map to the lineage-only label and can
    only earn lenient credit.
  - The top count was tied between 2+ labels in 38.7% (CellMarker) and 23.3% (PanglaoDB) of the 600
    scored Hao questions, usually k = 2; no overlap at all: 2.3% / 11.5% (all 1,200 questions).
- **CellTypist** was not used: it classifies individual cells from full expression profiles, not
  gene lists, so it can't be run on these questions or the difficulty knobs.

Results (Hao, same 600 questions per method, 300 per format, strict):

| | Ensembl IDs | Gene symbols |
|---|---|---|
| CellMarker lookup, tie-break | 28.7% | 28.7% |
| CellMarker lookup, fractional | 27.5% | 27.5% |
| Claude Sonnet 5 | 28.0% | 37.3% |
| Gemini 3.8 Flash | 22.0% | 37.7% |
| GPT-5.6 Terra | 21.3% | 32.7% |
| PanglaoDB lookup, tie-break | 16.0% | 16.0% |
| PanglaoDB lookup, fractional | 16.3% | 16.3% |
| Random guess | 3.3% | 3.3% |

(The lookups convert Ensembl IDs to symbols first, so their accuracy is the same in both formats.
Fractional numbers are over both formats together; the formats are identical for the lookups.)
Lenient: CellMarker 46.2% tie-break / 45.5% fractional; PanglaoDB 34.5% / 35.9%.
PBMC3k (8 types), strict: CellMarker 64.4% / 61.0%, PanglaoDB 36.9% / 41.9% (tie-break / fractional).

Statement of the gap: with gene symbols, the three LLMs score 4.0-9.0 points above the CellMarker
lookup (tie-break) and 5.2-10.2 points above it (fractional). With Ensembl IDs, Claude scores 0.7
points below the lookup (tie-break) and 0.5 points above it (fractional); Gemini and GPT score 5.5-7.3
points below it (5.5-6.2 fractional, 6.7-7.3 tie-break). Each LLM loses 9.3-15.7 points when the same genes are given as
Ensembl IDs instead of symbols; the lookup loses nothing, because it converts IDs to symbols first.

## Agent conditions (Hao, Claude only)

- **Subset** (`src/agent_subset.py`, seed 0, `data/agent_subset_hao.json`): 48 base questions, each in
  both formats (96 total, paired by format); 12 base questions per knob, spread evenly over levels;
  cell types dealt from a shuffled deck, so all 30 appear (2-4 questions each). Noise is included.
- **Agents** (`src/run_agent.py`): LangGraph prebuilt ReAct agent (`create_react_agent`, LangGraph
  1.2.12; marked deprecated in favour of `langchain.agents.create_agent` but functional) around
  `claude-sonnet-5` via the Anthropic API (langchain-anthropic 1.7.4), same prompt (`src/prompts.py`,
  shared with plain Claude), no temperature, default thinking, max_tokens 4000, same JSON parsing.
  At most 5 tool calls per question (later calls get a "budget used up" reply), 2-minute timeout.
  - generic: `web_search` (DuckDuckGo via langchain-community, top 5 results).
  - specialist: `gene_info` (MyGene.info: symbol/Ensembl -> symbol, name, summary);
    `cellmarker_gene_to_cell_types` and `cellmarker_cell_type_markers`, built from the same CellMarker
    marker lists as the lookup baseline (Hao 2021 excluded, top 50 per type, the 160 mapped names).
    The tools return CellMarker's own cell-type names, never the 30 Hao labels.
- Full trajectories: `results/agents/<condition>/<question id>.json`. Scored with `src/scoring.py`
  and the same judge (`src/score_agents.py`); one run per question.
- Tool quirk found: MyGene.info returns "TARP" for ENSG00000211689 (TRGC1, a gamma-delta T marker);
  in the one affected question the agent still answered correctly.

### Results on the 96-question subset (`score_agents.py`, one run each, n = 96: 48 symbol, 48 Ensembl)

| Condition | Strict | Lenient | Strict symbol | Strict Ensembl |
|---|---|---|---|---|
| plain Claude | 30.2% | 52.6% | 37.5% | 22.9% |
| CellMarker lookup | 27.1% | 46.9% | 27.1% | 27.1% |
| generic agent | 31.2% | 54.2% | 35.4% | 27.1% |
| specialist agent | 38.5% | 60.9% | 43.8% | 33.3% |
| specialist_nudge (post-hoc) | 42.7% | 64.6% | 41.7% | 43.8% |

- Paired (strict), specialist vs plain Claude: 10 vs 2 (sign test p = 0.039); vs CellMarker lookup:
  17 vs 6 (p = 0.035); vs generic agent: 7 vs 0 (p = 0.016).
- The generic agent called no tool on 84% of questions.
- The nudge row is post hoc (see below).

### Specialist on all 600 non-noise Hao questions

The original specialist (no nudge, unchanged) was then run on all 600 non-noise Hao questions - the set
plain Claude and the lookup were compared on - with `run_agent.py specialist --questions all_nonoise`
(the 72 non-noise subset questions already answered were reused, 528 new; 0 failures). Scored with
`score_agents.py --all-nonoise`; n = 600 per condition (300 symbol, 300 Ensembl), one run each.

| | Strict | Lenient | Strict symbol | Strict Ensembl |
|---|---|---|---|---|
| plain Claude | 32.7% | 54.7% | 37.3% | 28.0% |
| CellMarker lookup (tie-break) | 28.7% | 46.2% | 28.7% | 28.7% |
| specialist agent | 41.2% | 61.0% | 43.3% | 39.0% |

- Paired (strict): specialist right & plain Claude wrong on 58 questions, the reverse on 7 (two-sided
  sign test p = 4.3e-11); specialist right & lookup wrong on 90, the reverse on 15 (p = 3.3e-14).
- Strict by lineage (plain / lookup / specialist): B 44/32/50%, CD4 T 2/7/6%, CD8 T 8/0/12%, DC
  39/25/56%, Mono 78/50/85%, NK 38/27/43%, other 70/68/75%, other T 20/47/42%.
- Symbol vs Ensembl gap: specialist 4.3 points (43.3% vs 39.0%) vs 9.3 for plain Claude (37.3% vs 28.0%).
- 1.3 tool calls per question on average; $0.011 per question ($6.78 for the 600); 6.6 s mean latency.

### specialist_nudge — a SEPARATE condition, added AFTER the error analysis (post hoc)

- **Why it exists:** in the error analysis of the specialist's 59 wrong answers (96-question subset),
  21 were "no_cellmarker": the agent converted Ensembl IDs with `gene_info` and answered without
  querying CellMarker (20 of the 21 were Ensembl questions). The nudge tests whether telling it to
  query CellMarker fixes this. Because it was designed after seeing the specialist's errors on these
  same 96 questions, its results on them are post hoc and should not be read as an unbiased estimate.
- **Only difference from the specialist:** a one-line system prompt, "After converting any Ensembl IDs,
  always query CellMarker before answering." (The specialist and all other conditions have no system
  prompt.) Same tools, model, user prompt, limits and questions. The original specialist is unchanged.
- **Results** (same 96 questions, one run each, strict):

  | | Overall | Symbol | Ensembl | Lenient |
  |---|---|---|---|---|
  | specialist | 38.5% | 43.8% | 33.3% | 60.9% |
  | specialist_nudge | 42.7% | 41.7% | 43.8% | 64.6% |

  - Paired: nudge right & specialist wrong on 7 questions, the reverse on 3 (two-sided sign test
    p = 0.34 — not significant).
  - CellMarker was queried on 100% of questions with the nudge vs 40% (Ensembl) / 90% (symbol) without.
  - Of the specialist's 21 "no_cellmarker" errors, the nudge queried CellMarker on all 21 and answered
    **6 correctly** (strict); lenient mean on those 21 rose from 0.40 to 0.60.
  - Cost per question $0.016 vs $0.011; mean latency 8.0 s vs 5.9 s.

## Open decisions

- **"memory/effector":** your hand-check labelled "CD8+ memory/effector T cells (KLRB1+ GZMK+ subset,
  possibly MAIT-like)" as CD8 TEM (reading "memory/effector" as effector memory). Under the strict
  "X/Y alternatives are ignored" rule it would be CD8 T cell (subtype unclear), which is what the judge
  now answers. The test set keeps the hand-check label until this is decided.
- **Re-judge the "/" and "or" answers** under the current rule text (needs Anthropic credit, ~500
  judge calls) - see the "likely" rule above.

---

# Task 2: perturbation direction (Norman et al. 2019, CRISPRa in K562)

## Data
- scPerturb `NormanWeissman2019_filtered.h5ad` (Zenodo DOI 10.5281/zenodo.7041849, MD5 c870e696...,
  downloaded 2026-09-27): 111,445 K562 cells, 33,694 genes, raw counts. CRISPR **activation** screen.
  Only control (11,855 cells) and single-gene perturbations (105 genes, 57,831 cells) are used; the 131
  pair perturbations are not used yet. Three Norman labels use old symbols and are renamed to current
  HGNC symbols (C19orf26 -> CBARP, C3orf72 -> FOXL2NB, KIAA1804 -> MAP3K21; MAP3K21 is not measured).
- Normalisation: counts / cell total x 10,000, then natural log1p (cell totals equal the stored
  `ncounts`). Pseudobulk = mean log-normalised expression per perturbation and for control, computed
  by reading the gene-stored matrix in blocks of 2,000 genes (peak memory ~1.3 GB).
- CRISPRa check: the activated gene is higher than in control in 102/104 measured perturbations.

## Differential expression and labels (`src/task2_build.py`)
- **Fold change: log2FC = (mean log-normalised expression in X cells - mean in control cells) / ln 2**,
  i.e. the difference of mean log-expression, in log2 units. (This understates changes for lowly
  expressed genes compared with a log2 ratio of mean expression; e.g. the activated gene ATL1 is ~4x
  higher, log2 ratio 1.96, but difference-of-means log2FC 0.01.)
- **Test: Welch's t-test** on per-cell log-normalised values (X cells vs control cells), computed
  exactly from per-group sums and sums of squares; **Benjamini-Hochberg within each activated gene X**,
  across its tested target genes.
- Tested targets Y: detected (count > 0) in >= 10% of control cells (7,991 genes); Y != X.
- Labels: up = log2FC > 0.5 and padj < 0.05; down = log2FC < -0.5 and padj < 0.05;
  no_change = |log2FC| < 0.1 and padj > 0.5; everything else dropped. Over all 839,001 tested pairs:
  932 up, 1,299 down, 414,483 no_change, 422,287 dropped.

## Question set (`data/task2_questions.json`, 200 questions, seed 0)
- 67 up, 67 down, 66 no_change; all 105 activated genes (max 2 questions each; cap 5); 159 target genes
  (cap 3 per target); targets with clone-style names (RP11-..., AC012345.1, ...) excluded.
- **Expression matching:** up and down are sampled first; each no_change question is matched to a random
  up/down question on the target's control detection rate (within 0.02, else 0.05; all 66 matched).
  Without matching, no_change targets had a median detection of 32% vs 94-96% for up/down. After
  matching, no_change vs up+down detection: KS D = 0.098, p = 0.74.
- **Remaining imbalance (not corrected):** up targets are less expressed than down targets (median
  detection 75% vs 98%; KS p < 0.001); detection rate alone separates down from up with AUC 0.80. This
  reflects the biology (activation switches on quieter genes and represses highly expressed ones) but
  is a possible shortcut for the direction questions.
- Question text: "In K562 cells (a human leukemia cell line), gene X is activated with CRISPR activation
  (CRISPRa). What happens to the expression of gene Y? Answer up, down, or no_change, as JSON with a
  confidence 0-1." Gene symbols only.

## Baselines (`src/task2_baselines.py`, no API calls)
- always no_change.
- co-expression sign: Pearson r of X and Y across control cells only; r > 0 -> up, r < 0 -> down;
  X undetected in controls -> no_change (4 questions).
- co-expression with |r| <= 0.05 -> no_change: **threshold set when the script was first written,
  before any result was seen; not tuned.**
- CollecTRI (via decoupler 2.2.0 / OmniPath, snapshot `data/task2/collectri_human.csv`, 42,990 links,
  downloaded 2026-09-27): X -> Y link weight +1 -> up, -1 -> down; no link -> no_change. Note that
  CollecTRI labels links without sign evidence as "default activation" (+1).
- Metrics: 3-class accuracy, macro-F1, and direction accuracy = share of up/down questions predicted
  with the correct direction (no_change predictions count as wrong).

| Baseline | Accuracy | Macro-F1 | Direction acc. (up/down) |
|---|---|---|---|
| always no_change | 33.0% | 0.165 | 0.0% |
| co-expression sign | 42.0% | 0.345 | 61.9% |
| co-expression, \|r\| <= 0.05 -> no_change (fixed in advance) | 39.0% | 0.298 | 11.2% |
| CollecTRI | 33.5% | 0.176 | 0.7% |

CollecTRI coverage: 1/200 questions (0.5%) have an X -> Y link (81/200 have an X with any CollecTRI
targets).

## Task 2 specialist agent
`src/run_agent_task2.py`: same setup as the Task 1 agents (claude-sonnet-5, Anthropic API, LangGraph
prebuilt ReAct agent, question text as the only prompt, no system prompt, max 5 tool calls, 2-minute
timeout, full trajectories saved). Three tools, none of which can see perturbation outcomes:
1. `gene_info` - MyGene.info (symbol/Ensembl -> symbol, name, summary).
2. `coexpression_in_control` - Pearson r of X and Y across the 11,855 control cells (from
   `data/norman2019/task2_ctrl_corr.npz`, written by `task2_build.py`), plus the share of control
   cells expressing each gene. Returns exactly the numbers the co-expression baseline uses (checked on
   all 200 questions).
3. `collectri_lookup` - the same CollecTRI snapshot as the baseline: whether X is a TF in CollecTRI,
   its link to Y (sign, "default activation" flag, number of references), up to 30 other targets.
Answers are parsed as the first up / down / no_change value in the JSON reply.

- Run on all 200 questions after a 5-question pilot (0 failures, $1.55, 2.1 tool calls per question on
  average; co-expression queried for (X, Y) on 98% of questions, CollecTRI on 77%, gene_info on 32%).

## Task 2 results (n = 200, one run per question)

Plain models: `run_eval.py --dataset task2` - Claude Sonnet 5 via the Anthropic API; GPT-5.6 Terra and
Gemini 3.8 Flash via OpenRouter (Gemini's Google API billing was unavailable); the question text as the
only prompt; same parser for every method (`src/task2_common.py`). Gemini's one API error (t2_004) was
retried once and answered. Scored with `src/task2_score.py`; unparsable answers count as wrong.
Direction accuracy is on the 134 up/down questions; matched pairs are the 30 expression-matched pairs
(60 questions) fixed before any model run.

| Method | Accuracy | Macro-F1 | Direction, all up/down | Direction, matched pairs | no_change on up/down |
|---|---|---|---|---|---|
| baseline: always no_change | 33.0% | 0.165 | 0.0% | 0.0% | 100% |
| baseline: co-expression sign | 42.0% | 0.345 | 61.9% | 60.0% | 2% |
| baseline: co-expression, \|r\| <= 0.05 -> no_change (fixed in advance) | 39.0% | 0.298 | 11.2% | 10.0% | 88% |
| baseline: CollecTRI | 33.5% | 0.176 | 0.7% | 0.0% | 99% |
| baseline: expression level only (threshold fixed in advance) | 49.0% | 0.391 | 73.1% | 50.0% | 0% |
| Claude Sonnet 5 | 41.5% | 0.359 | 20.9% | 15.0% | 64% |
| GPT-5.6 Terra | 46.0% | 0.430 | 29.9% | 31.7% | 49% |
| Gemini 3.8 Flash | 45.5% | 0.410 | 26.1% | 25.0% | 58% |
| specialist agent (Claude + 3 tools) | 39.0% | 0.288 | 9.7% | 1.7% | 88% |

- **Agent vs plain Claude (paired, same questions):** 3-class, agent right & plain wrong 14 vs the
  reverse 19 (two-sided sign test p = 0.49); direction on up/down, 4 vs 19 (p = 0.003). The agent gave
  the same answer as plain Claude on 75% of questions and answered no_change 183/200 times (plain: 141).
- **Agent wrong answers (122/200), classified from the co-expression value in its own tool output**
  ("weak" = |r| <= 0.05, the baseline's threshold fixed before any model run):
  weak co-expression read as no_change 105; |r| > 0.05 but still no_change 10; no_change without
  querying co-expression 3; wrong direction following the sign of r 2; wrong direction against the sign
  1; claimed an effect on a no_change question 1. Control-cell co-expression was weak (|r| <= 0.05) for
  90% of queried pairs whatever the truth (53 up, 62 down, 62 no_change), and the agent answered
  no_change on 94% of those. When |r| > 0.05 and it did answer up/down, it always followed the sign of r.
- **Extra metric, defined post hoc (after seeing the model results): direction accuracy when committed**
  = among up/down questions where the method answered up or down, the share with the correct direction.
  It separates "declined to commit" (no_change) from "chose the wrong direction". Not a pre-planned
  metric; small n for methods that rarely commit.

  | Method | All up/down: committed, correct | Matched pairs: committed, correct |
  |---|---|---|
  | Claude Sonnet 5 | 48/134, 58.3% | 21/60, 42.9% |
  | GPT-5.6 Terra | 69/134, 58.0% | 33/60, 57.6% |
  | Gemini 3.8 Flash | 57/134, 61.4% | 26/60, 57.7% |
  | specialist agent | 16/134, 81.2% | 2/60, 50.0% |
  | baseline: co-expression sign | 131/134, 63.4% | 59/60, 61.0% |
  | baseline: expression level only | 134/134, 73.1% | 60/60, 50.0% |
- Direction biases (what the models answered on true up / true down questions): Claude answered up
  21 / 17 times and down 3 / 7; GPT and Gemini answered down more than up in both.

---

# Cost (logged per answer in the results files, US$)

| Component | Cost |
|---|---|
| Task 1, PBMC3k pilot (3 models x 3 runs) | $10.74 |
| Task 1, Hao, plain models | $8.65 |
| Task 1, agents (generic, specialist incl. all 600, nudge) | $9.19 |
| Task 2, plain models | $1.80 |
| Task 2, specialist agent | $1.55 |
| **Total logged** | **$31.94** |

Not logged: LLM-judge calls, one-question tests, and runs discarded before cost logging was added.

**Account totals (actual spend): [PLACEHOLDER - to be filled in from the Anthropic and OpenRouter accounts]**
