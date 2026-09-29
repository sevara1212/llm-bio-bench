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

- **Judge route (provenance):** from 2026-09-28 the Anthropic API credit was exhausted, and new judge
  decisions were made with the SAME judge (`claude-sonnet-5`, same instructions, label list and JSON
  schema) called through OpenRouter (`anthropic/claude-sonnet-5`, structured output, reasoning effort
  low), paced to ~15 requests/minute. On the 41-case test set this route scored 40/41, identical to the
  Anthropic route (the miss is the open "memory/effector" case). 176 of the 1,290 cached decisions were
  made this way (all of them for answers first seen in these later runs; no earlier decision was re-made);
  each is marked in `results/judge_cache_hao_provenance.json`, the others came from the Anthropic API.

## Judge check: second pass over the judge cache (`src/judge_check.py`, run after all results)
- Every cached judge decision was judged a second time, independently: same judge model (claude-sonnet-5),
  instructions, label list and parsing, called through OpenRouter. New labels are in
  `results/judge_check_hao.json`; the main cache and all reported scores still use the first pass.
  Stopped at the $4 OpenRouter balance floor after **1,140 of 1,290 decisions** (981 of the 1,114 first made
  via the Anthropic API, 159 of the 176 made via OpenRouter; no failed calls; the balance fell by $5.63).
- **Agreement:** exact label 90.9%, same lineage 97.4%. Decisions first made via OpenRouter (recently, under the
  final rule text): 99.4% exact, 100% lineage. Decisions first made via the Anthropic API: 89.5% exact, 96.9%
  lineage (88.2% for answers with "/" or "or" alternatives, 91.8% without). The two routes are confounded with
  time: many Anthropic-route labels were made before later rule changes (see "Rule changes" above; the 509
  alternative-wording answers were never re-judged), so this gap cannot be attributed to the route alone.
- **Disagreements** are mostly ambiguous memory/effector wordings that the second pass maps to "subtype
  unclear": CD8 TEM -> CD8 T cell (subtype unclear) 28, CD4 TCM -> CD4 T cell (subtype unclear) 15, CD8 TCM
  -> CD8 T cell (subtype unclear) 11 (list in `results/judge_check_disagreements.csv`).
- **Effect on the headline Task 1 scores** (600 questions, strict, first -> second pass): Claude 32.7% ->
  31.2%, Gemini 29.8% -> 30.0%, CellMarker 28.7% (rule-matched, unchanged), GPT 27.0% -> 26.8%, PanglaoDB
  16.0% (unchanged); specialist agent 41.2% -> 39.8%. Specialist vs plain Claude, paired: 58 vs 7 ->
  58 vs 6 (p = 9.0e-12). The ranking and every conclusion are unchanged.
- **Effect on Task 3 on Hao:** mean strict accuracy of runs with labels 40.5% -> 40.4% (15 steps), 47.7% ->
  48.0% (30 steps), 25.5% -> 25.6% (stateful); the largest change in any run is 2.9 points.

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

### Specialist agent with GPT-5.6 Terra and Gemini 3.8 Flash (96-question subset)

- Same agent as the Claude specialist (tools, prompt, 5-call budget, 2-minute timeout), with GPT or
  Gemini through OpenRouter at temperature 0 (as in their plain runs) and 2 parallel workers (new
  OpenRouter accounts are limited to 20 requests/minute per model; rate-limited questions were retried).
  `run_agent.py specialist --model openai/gpt-5.6-terra` / `google/gemini-3.8-flash`.
- Plain GPT and Gemini had not answered the subset's noise questions (they were run with noise skipped);
  those 23 / 22 questions were added with `run_eval.py ... --ids data/agent_subset_hao.json` ($0.24), so
  every comparison uses all 96 questions.

| Condition (n = 96) | Strict | Lenient | Strict symbol | Strict Ensembl | Tool calls / question | Cost / question |
|---|---|---|---|---|---|---|
| plain GPT-5.6 Terra | 29.2% | 51.0% | 37.5% | 20.8% | - | - |
| specialist agent (GPT) | 43.8% | 63.0% | 41.7% | 45.8% | 1.22 (31% none) | $0.005 |
| plain Gemini 3.8 Flash | 30.2% | 53.1% | 35.4% | 25.0% | - | - |
| specialist agent (Gemini) | 38.5% | 59.9% | 39.6% | 37.5% | 2.75 (11.5% hit the budget) | $0.005 |
| specialist agent (Claude), for reference | 38.5% | 60.9% | 43.8% | 33.3% | 1.24 | $0.011 |

- Paired (strict) vs each model's own plain answers: GPT agent right & plain wrong on 15 questions, the
  reverse on 1 (sign test p = 0.0005); Gemini 11 vs 3 (p = 0.057). Vs Claude's agent: GPT 8 vs 3
  (p = 0.227); Gemini 4 vs 4 (p = 1.000).

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
| Gemini 3.8 Flash | 45.5% | 0.410 | 26.1% | 25.0% | 57% |
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

## Post hoc condition A (SEPARATE, added after the Claude agent results): the Task 2 specialist agent with GPT and Gemini
- Same tools, question text, limits (5 tool calls, 2-minute timeout) and parser as the Claude agent; models
  GPT-5.6 Terra and Gemini 3.8 Flash via OpenRouter at temperature 0 (`run_agent_task2.py --model openai/...`
  / `google/...`), all 200 questions, one run each. 3 GPT questions hit OpenRouter rate limits (429) and were
  rerun; every question finished. Gemini gave 4 replies with no parsable answer (counted as wrong).
  Analysis: `task2_agent_analysis.py --model <model>`.

| Method (post hoc rows marked) | Accuracy | Macro-F1 | Direction, all up/down | Direction, matched pairs | no_change on up/down |
|---|---|---|---|---|---|
| GPT-5.6 Terra, plain | 46.0% | 0.430 | 29.9% | 31.7% | 49% |
| specialist agent, GPT + 3 tools (post hoc) | 43.5% | 0.370 | 18.7% | 20.0% | 72% |
| Gemini 3.8 Flash, plain | 45.5% | 0.410 | 26.1% | 25.0% | 57% |
| specialist agent, Gemini + 3 tools (post hoc) | 39.5% | 0.312 | 11.9% | 6.7% | 81% |
| specialist agent, Claude + 3 tools (main) | 39.0% | 0.288 | 9.7% | 1.7% | 88% |

- **Vs its own plain answers (paired, sign test):** GPT agent 3-class 17 right where plain was wrong vs 22
  the reverse (p = 0.522); direction 6 vs 21 (p = 0.006). Gemini agent 3-class 11 vs 23 (p = 0.058);
  direction 2 vs 21 (p = 6.6e-05). With tools both models committed to a direction less often and lost
  direction accuracy, as Claude did (4 vs 19, p = 0.003).
- **Vs Claude's agent:** GPT agent 3-class 16 vs 7 (p = 0.093), direction 16 vs 4 (p = 0.012); Gemini agent
  3-class 6 vs 5 (p = 1.0), direction 6 vs 3 (p = 0.508).
- **Weak co-expression read as no_change - yes, for both.** When the (X, Y) control-cell correlation it
  queried was weak (|r| <= 0.05), the GPT agent answered no_change 85% of the time (n = 162) and the Gemini
  agent 89% (n = 177; Claude's agent 94%). This was the largest class of wrong answers: 81 of GPT's 113 and
  98 of Gemini's 121 (Claude's agent: 105 of 122). When |r| > 0.05 and they answered up/down, they
  followed the sign of r 89% (GPT) and 100% (Gemini) of the time.
- Tool use: mean 2.17 (GPT) and 2.83 (Gemini) tool calls per question; co-expression queried for (X, Y) on
  90% / 98% of questions, CollecTRI on 68% / 75%, gene_info on 57% / 90%.

---

# Task 3: end-to-end analysis by an agent (PBMC3k, raw counts)

## Setup (`src/task3_agent.py`, `src/task3_sandbox.py`)
- **Input:** each run gets a fresh folder containing only the raw PBMC3k counts as a 10x
  `filtered_gene_bc_matrices/hg19/` folder (matrix.mtx, genes.tsv, barcodes.tsv; 2,700 cells x 32,738 genes),
  written from scanpy's cached `pbmc3k_raw.h5ad` - the same raw input `prepare_data.py` starts from.
- **Agent:** LangGraph prebuilt ReAct agent (`create_react_agent`) with one tool, `run_python(code) -> output`.
  All three models via OpenRouter (`anthropic/claude-sonnet-5`, `openai/gpt-5.6-terra`,
  `google/gemini-3.8-flash`; provider and model recorded per run), max_tokens 8,000, no temperature set,
  no system prompt. The only instruction, identical for every model: "This folder contains a raw single-cell RNA-seq count matrix from human peripheral blood mononuclear cells (PBMCs). Perform standard quality control, normalisation and clustering, identify marker genes, and assign a cell type label to every cluster. Save labels.csv with columns barcode,cell_type. Work step by step and check your outputs."
- **Tool description** (also identical): each call runs in a new Python process with a 60-second limit,
  variables do not persist (save intermediate results to files), scanpy / anndata / pandas / numpy / scipy /
  leidenalg installed, no internet. Each output ends with "[execution k/15, exit code X; n left]".
- **Limits:** at most 15 code executions ("steps") per run; a run stops at $1.50 (checked after every model
  reply, from token counts at OpenRouter list prices); before every run the live OpenRouter balance was
  checked and all runs would stop below $3 (never triggered). A run without labels.csv counts as failed.
- **Runs:** 1 pilot run per model, then 6 more (7 per model, 21 in total), interleaved across models, each
  saved as soon as it ended.

## Sandbox rules
- macOS `sandbox-exec` (Seatbelt) profile, deny by default: read access only to the Python installation,
  `/opt/homebrew` (installed software and its libraries), the project's `venv` and the run's own folder;
  write access only inside the run folder; all network access denied; a fresh process per execution with a
  minimal environment (no API keys, HOME = the run folder); the process group is killed after 60 s.
- Tested before any model run: reading the answer key, the project `.env`, the home and project folders,
  other runs' folders, running `/bin/cat` on the answer key, writing outside the folder, raw sockets, HTTPS
  and DNS were all blocked; the scanpy pipeline ran (15 s).
- Run folders are under `/private/tmp/llmbio_task3_runs/`, outside the project: anndata's settings loader
  searches parent folders for a `.env` file and, inside the project, found the API-key file (the read was
  blocked, but it crashed the import).
- Output shown to the agent is truncated to ~50 lines / 4,000 characters (first 20 + last 30 lines); the
  untruncated output (up to 20,000 characters) is kept in the trajectory. joblib prints a one-line warning
  in every execution because the sandbox blocks its shared-memory semaphores (it falls back to serial mode).
- Trajectories: `results/task3/runs/<run_id>/trajectory.json` (every code cell, its full and shown output,
  the model's text, tokens, cost), plus `code/`, `labels.csv` and `agent_obs_*.csv` (the obs table of every
  .h5ad the agent saved). Hidden reasoning was not returned by OpenRouter for these models.

## Scoring (`src/task3_score.py`)
- **Answer key:** the per-cell labels of this project's PBMC3k pipeline (`prepare_data.py`, saved as
  `data/pbmc3k_expert_cells.csv`: 2,638 cells after the tutorial QC, 8 Leiden clusters with tutorial cell
  types). These are tutorial-derived labels, not an independent expert annotation.
- Only cells present in both the agent's labels.csv and the answer key are scored; coverage = share of the
  2,638 expert cells the agent kept (100% in every successful run).
- Free-text labels are mapped to the 8 types with the Task 1 PBMC3k keyword rules (`normalize_pbmc`,
  `score_pbmc`): strict = exact type, lenient = 0.5 for "T cell" on CD4/CD8 or "monocyte" on CD14/FCGR3A.
- **ARI:** adjusted Rand index of the partition given by the agent's labels vs the 8 expert Leiden clusters
  (splitting a type, e.g. CD4 naive vs memory, lowers ARI even when both labels are correct).
- **Diagnosis:** for every expert type < 80% correct, the wrongly labelled cells are located in the
  agent's own clusters (the saved cluster column that best matches its labels.csv): wrong cells in clusters
  dominated by another type, or in a mixed cluster (second type >= 25%), = merged in clustering; wrong cells
  in a pure cluster of their own type = right cluster, wrong name. Runs that saved no clusters, and the
  failed run, were classified by reading their trajectories.

## Results (21 runs)

| Run | Status | Cells kept | Strict | Lenient | ARI | Labels | Steps | Code errors | Cost (est.) | Time | Where it went wrong |
|---|---|---|---|---|---|---|---|---|---|---|---|
| claude_run01 | ok | 2638 | 97.1% | 97.1% | 0.94 | 8 | 14 | 2 | $0.225 | 139 s | none (all types >= 80% correct) |
| gpt_run01 | ok | 2700 | 88.7% | 88.7% | 0.58 | 9 | 11 | 1 | $0.186 | 116 s | clustering: switched to Leiden resolution 1.0 at step 6 (at 0.6 the DCs were separate), merging all 36 DCs into the non-classical monocyte cluster |
| gemini_run01 | ok | 2638 | 93.1% | 93.1% | 0.63 | 9 | 15 | 1 | $0.080 | 140 s | clustering: 85/339 CD8 T cells in a CD4-majority cluster |
| claude_run02 | ok | 2638 | 93.0% | 93.0% | 0.87 | 8 | 13 | 0 | $0.231 | 178 s | clustering: 71/339 CD8 T cells in a CD4-majority cluster |
| gpt_run02 | ok | 2700 | 93.5% | 93.5% | 0.87 | 8 | 8 | 0 | $0.129 | 124 s | clustering: 77/339 CD8 T cells in a CD4-majority cluster |
| gemini_run02 | ok | 2638 | 93.0% | 93.0% | 0.87 | 8 | 15 | 4 | $0.063 | 125 s | clustering: 71/339 CD8 T cells in a CD4-majority cluster |
| claude_run03 | ok | 2638 | 97.1% | 97.1% | 0.94 | 8 | 10 | 0 | $0.138 | 120 s | none |
| gpt_run03 | ok | 2700 | 91.6% | 91.6% | 0.66 | 8 | 7 | 0 | $0.093 | 126 s | clustering: 110/339 CD8 T cells in a CD4-majority cluster; all 13 megakaryocytes merged into the FCGR3A monocyte cluster |
| gemini_run03 | ok | 2638 | 20.9% | 55.5% | 0.82 | 6 | 15 | 3 | $0.078 | 187 s | step budget: re-ran the pipeline in most steps (3 NameErrors, state assumed to persist); final coarse clustering (resolution 0.5, 6 clusters) at step 14 merged CD8 T + NK ('NK cells') and left CD4 as generic 'T cells', monocytes unsplit |
| claude_run04 | ok | 2638 | 97.1% | 97.1% | 0.94 | 8 | 10 | 0 | $0.138 | 118 s | none |
| gpt_run04 | ok | 2700 | 94.2% | 94.2% | 0.65 | 10 | 11 | 1 | $0.187 | 144 s | clustering: 87/339 CD8 T cells in a CD4-majority cluster |
| gemini_run04 | ok | 2638 | 42.2% | 64.5% | 0.84 | 7 | 15 | 4 | $0.080 | 184 s | step budget: 3 NameErrors, labels written at step 15; CD8 T + NK merged ('NK cells'), CD4 left as generic 'T cells' |
| claude_run05 | ok | 2638 | 97.1% | 97.1% | 0.94 | 8 | 11 | 1 | $0.154 | 116 s | none |
| gpt_run05 | ok | 2700 | 92.6% | 92.6% | 0.85 | 8 | 11 | 0 | $0.173 | 105 s | clustering: 85/339 CD8 T cells in a CD4-majority cluster |
| gemini_run05 | ok | 2638 | 84.2% | 84.2% | 0.84 | 7 | 15 | 5 | $0.066 | 127 s | clustering: CD8 T + NK merged in one cluster (59% / 37%), named 'NK cells'; labels written at step 15 |
| claude_run06 | ok | 2638 | 93.0% | 93.0% | 0.87 | 8 | 14 | 0 | $0.245 | 156 s | clustering: 71/339 CD8 T cells in a CD4-majority cluster |
| gpt_run06 | ok | 2700 | 84.3% | 84.3% | 0.84 | 7 | 10 | 2 | $0.125 | 99 s | clustering: CD8 T + NK merged in one cluster (58% / 35%), named 'NK cell' from shared cytotoxic markers |
| gemini_run06 | failed | - | 0.0% | 0.0% | - | - | 15 | 4 | $0.098 | 164 s | ran out of steps: 15/15 used, 4 code errors (incl. NameErrors), no labels.csv; its final message nevertheless claimed the analysis was complete |
| claude_run07 | ok | 2638 | 92.5% | 92.5% | 0.62 | 9 | 13 | 0 | $0.208 | 138 s | clustering: 84/339 CD8 T cells in a CD4-majority cluster |
| gpt_run07 | ok | 2700 | 92.9% | 92.9% | 0.78 | 9 | 10 | 2 | $0.143 | 90 s | clustering: 77/339 CD8 T cells in a CD4-majority cluster |
| gemini_run07 | ok | 2638 | 44.0% | 66.2% | 0.86 | 7 | 15 | 5 | $0.073 | 150 s | step budget: 3 NameErrors, labels written at step 15; CD8 T + NK merged ('NK cells'), CD4 left as generic 'T cells' |

### Consistency across runs (7 per model; the failed run counts as 0)

| Model | Runs with labels | Strict mean (SD) | Strict range | Lenient mean | ARI mean (SD) | ARI range | Labels range | Steps mean | Code errors per run | Cost per run (est.) |\n|---|---|---|---|---|---|---|---|---|---|---|\n| Claude Sonnet 5 | 7/7 | 95.3% (2.3%) | 92.5%-97.1% | 95.3% | 0.87 (0.12) | 0.62-0.94 | 8-9 | 12.1 | 0.4 | $0.191 |\n| GPT-5.6 Terra | 7/7 | 91.1% (3.5%) | 84.3%-94.2% | 91.1% | 0.75 (0.12) | 0.58-0.87 | 7-10 | 9.7 | 0.9 | $0.148 |\n| Gemini 3.8 Flash | 6/7 | 53.9% (37.0%) | 0.0%-93.1% | 65.2% | 0.81 (0.09) | 0.63-0.87 | 6-9 | 15.0 | 3.7 | $0.077 |

## Failure classes
- **Clustering - CD8 T cells partly merged into a CD4-majority cluster:** 10 of 21 runs (Claude 3, GPT 5,
  Gemini 2), 71-110 of the 339 CD8 T cells each time.
- **Clustering - CD8 T and NK cells merged into one cluster, named NK:** 5 runs (GPT 1, Gemini 4). The
  mixed cluster (~430-445 cells, about 58-59% CD8 T / 35-37% NK where the agent's clusters were saved) was
  named after the cytotoxic markers the two types share (NKG7, GZMA, CST7, CCL5), so ~250-270 CD8 T cells
  were called NK.
- **Clustering - other merges:** dendritic cells into non-classical monocytes (GPT run 1, caused by the
  choice of resolution 1.0 at step 6); megakaryocytes into FCGR3A monocytes (GPT run 3).
- **Labelling too coarse:** 3 Gemini runs labelled all CD4 T cells as generic "T cells" (lenient credit
  only), in each case in a final pipeline written at step 14-15.
- **Ran out of steps / did not keep state:** Gemini used all 15 steps in every run and had 26 code errors
  in total (3.7 per run, mostly `NameError`s from assuming variables persisted between calls, despite the
  tool description), re-running the whole pipeline in most steps. One run (Gemini run 6) produced no
  labels.csv and still ended with a message claiming the analysis was complete.
- **QC:** never too strict (coverage 100%). GPT kept all 2,700 cells in every run (62 more than the
  tutorial QC; e.g. run 1 removed 2 cells and re-assigned them by nearest neighbours) - looser than the
  answer key but with no effect on accuracy here, because only shared cells are scored. Claude and Gemini
  kept exactly the tutorial's 2,638 cells in every successful run (the tutorial thresholds: 200-2,500 genes,
  < 5% mitochondrial counts).
- **Normalisation and marker finding:** no failures traced to these. Every run used total-count
  normalisation to 10,000 followed by log1p; where types were merged or misnamed, the markers the agent
  printed were consistent with the (merged) cluster it had made.

## Cost
- Estimated from token counts at OpenRouter list prices: $2.91 for the 21 runs (Claude $1.34, GPT $1.04,
  Gemini $0.54).
- Actual OpenRouter balance change over the pilot and the 18 runs: $28.74 -> $26.41 ($2.33); the estimate
  ignores provider prompt-caching discounts.

## Limitations
- **PBMC3k is a famous tutorial dataset** (the scanpy and Seurat introductory tutorial) and the answer key
  is itself the tutorial pipeline's output; models may have memorised the workflow and its cell types.
  Claude and Gemini reproduced the tutorial QC thresholds exactly in every successful run.
- One dataset, one instruction, 7 runs per model; no temperature control; the step and cost limits shape
  behaviour (Gemini's failures are mostly step-budget failures).
- ARI compares with 8 tutorial clusters, so reasonable finer splits are penalised.

---

# Task 3 on Hao 2021 (30 fine cell types)

## Setup (`src/task3_hao_data.py`, `task3_agent.py --dataset hao`, `src/task3_hao_score.py`)
- **Agent's input:** `raw_counts.h5ad` - raw RNA counts of up to 300 cells per Hao celltype.l2 type (30 types,
  doublets dropped, seed 0; 7,841 cells x 20,264 genes; 76-300 cells per type), gene symbols as var_names
  plus Ensembl gene_ids, cell barcodes. No labels, no protein (ADT) data, no embeddings (the CELLxGENE file's
  `X_apca`/`X_spca`/WNN embeddings are not copied), no other metadata; cells shuffled (seed 0). Counts
  verified identical to Hao's `raw.X`. It is an .h5ad rather than a 10x .mtx folder so that loading takes
  0.5 s instead of a large share of the 60-second execution limit.
- **Same** sandbox, instruction, tool description (fresh process per call), 15-step cap, $1.50 per run and
  models (all via OpenRouter) as PBMC3k Task 3; 7 runs per model. The balance floor was $4-5 and never
  reached.
- **Answer key** `data/task3_hao_expert_cells.csv` (barcode, celltype.l1, celltype.l2), outside the sandbox.
- **Scoring:** agent labels -> Hao labels with the Task 1 rules (`judge.rule_match`, then the LLM judge; see
  "Judge route (provenance)"). strict_fine = exact celltype.l2; lenient_fine = Task 1 lenient score; coarse =
  predicted lineage equals celltype.l1 (fine labels by their l1, lineage-only labels by their lineage, "T cell
  (unclear)" counts for CD4 T, CD8 T and other T); ARI of the agent's label partition vs celltype.l2 / l1.
- **"Proliferating lymphocyte (lineage unclear)" counts as wrong at the lineage level** (it names no lineage).
  It affects 5,852 of 125,337 scored cells (4.7%) in the 16 successful 15-step runs - 5,734 of them truly
  proliferating cells (NK Proliferating 3,244, CD4 Proliferating 1,352, CD8 Proliferating 1,138); 4,209 of
  78,398 (5.4%) in the 30-step runs and 2,962 of 47,046 (6.3%) in the stateful runs.
- **Tutorial thresholds:** detected from each run's code (max genes 2,500 AND mitochondrial < 5%).

## Results (21 runs, 15 steps)

| Run | Status | Cells kept | Strict (30 types) | Lenient | Coarse (lineage) | ARI (l2) | Labels | Steps | Errors | Cost (est.) | Stop reason |
|---|---|---|---|---|---|---|---|---|---|---|---|
| hao_claude_run01 | ok | 7771 | 53.1% | 66.9% | 74.6% | 0.58 | 22 | 15 | 0 | $0.493 | agent finished (all steps used) |
| hao_gpt_run01 | ok | 7841 | 42.0% | 58.4% | 74.8% | 0.51 | 20 | 10 | 1 | $0.138 | agent finished |
| hao_gemini_run01 | ok | 7841 | 32.6% | 54.1% | 69.2% | 0.48 | 16 | 15 | 3 | $0.084 | agent finished (all steps used) |
| hao_claude_run02 | failed | - | 0.0% | 0.0% | 0.0% | - | - | 15 | 0 | $0.370 | step limit reached |
| hao_gpt_run02 | ok | 7841 | 46.6% | 61.6% | 70.4% | 0.54 | 21 | 10 | 1 | $0.185 | agent finished |
| hao_gemini_run02 | ok | 7841 | 29.0% | 49.1% | 63.1% | 0.44 | 16 | 15 | 2 | $0.096 | agent finished (all steps used) |
| hao_claude_run03 | ok | 7838 | 45.6% | 59.1% | 66.4% | 0.53 | 20 | 15 | 0 | $0.347 | agent finished (all steps used) |
| hao_gpt_run03 | ok | 7841 | 37.5% | 56.8% | 76.2% | 0.51 | 19 | 13 | 2 | $0.200 | agent finished |
| hao_gemini_run03 | ok | 7841 | 46.5% | 61.9% | 71.0% | 0.53 | 19 | 15 | 3 | $0.096 | agent finished (all steps used) |
| hao_claude_run04 | ok | 7795 | 50.4% | 63.5% | 71.7% | 0.50 | 20 | 15 | 0 | $0.458 | agent finished (all steps used) |
| hao_gpt_run04 | ok | 7841 | 40.2% | 58.0% | 76.8% | 0.51 | 20 | 8 | 1 | $0.122 | agent finished |
| hao_gemini_run04 | ok | 7841 | 35.7% | 55.6% | 69.1% | 0.48 | 16 | 15 | 3 | $0.079 | agent finished (all steps used) |
| hao_claude_run05 | failed | - | 0.0% | 0.0% | 0.0% | - | - | 15 | 0 | $0.417 | step limit reached |
| hao_gpt_run05 | ok | 7841 | 39.1% | 56.6% | 74.3% | 0.50 | 20 | 10 | 1 | $0.143 | agent finished |
| hao_gemini_run05 | ok | 7841 | 25.2% | 48.6% | 65.7% | 0.43 | 15 | 15 | 3 | $0.087 | agent finished (all steps used) |
| hao_claude_run06 | failed | - | 0.0% | 0.0% | 0.0% | - | - | 15 | 0 | $0.298 | step limit reached |
| hao_gpt_run06 | ok | 7841 | 42.9% | 55.8% | 66.2% | 0.53 | 22 | 12 | 1 | $0.242 | agent finished |
| hao_gemini_run06 | ok | 7841 | 37.8% | 59.0% | 77.6% | 0.46 | 19 | 15 | 1 | $0.119 | agent finished (all steps used) |
| hao_claude_run07 | failed | - | 0.0% | 0.0% | 0.0% | - | - | 15 | 1 | $0.349 | step limit reached |
| hao_gpt_run07 | ok | 7841 | 43.7% | 60.6% | 75.2% | 0.60 | 22 | 11 | 0 | $0.230 | agent finished |
| hao_gemini_run07 | failed | - | 0.0% | 0.0% | 0.0% | - | - | 15 | 3 | $0.073 | agent finished (all steps used) |

| Model | Runs with labels | Strict mean, all runs (SD) | Strict, runs with labels: mean (range) | Coarse, runs with labels | ARI (l2), runs with labels | Steps mean (range) | Errors / run | Cost / run (est.) |
|---|---|---|---|---|---|---|---|---|
| Claude Sonnet 5 | 3/7 | 21.3% (26.7%) | 49.7% (45.6%-53.1%) | 70.9% | 0.54 | 15.0 (15-15) | 0.1 | $0.390 |
| GPT-5.6 Terra | 7/7 | 41.7% (3.1%) | 41.7% (37.5%-46.6%) | 73.4% | 0.53 | 10.6 (8-13) | 1.0 | $0.180 |
| Gemini 3.8 Flash | 6/7 | 29.5% (14.7%) | 34.5% (25.2%-46.5%) | 69.3% | 0.47 | 15.0 (15-15) | 2.6 | $0.090 |

- **Claude hit the 15-step limit without saving labels.csv in 4 of 7 runs**; its 3 finished runs are the most
  accurate (45.6-53.1% strict). GPT finished every run in 8-13 steps. Gemini used all 15 steps in every run.
- **Main error, all models:** T-cell subtypes. More than half of CD8 T cells were given CD4-lineage labels (15 steps:
  Claude 55%, GPT 60%, Gemini 62% of CD8 T cells) and most "other T" cells (MAIT, gdT, dnT) were called CD4
  or CD8 (Claude 76%, GPT 55%, Gemini 70%); CD4 T subtypes were exactly right for only 14-19% of cells. B cells, DCs and
  monocytes were almost always right at the lineage level.
- **Tutorial thresholds were not reused on Hao** (0/21 runs): typical choices were < 5,500-6,000 genes and < 15% mitochondrial (Claude)
  and >= 200-500 genes with < 12-15% mitochondrial (GPT). Coverage was 99-100%.
- **Compared with PBMC3k Task 3** (8 types): strict accuracy of runs with labels fell from 95.3% to 49.7%
  (Claude), 91.1% to 41.7% (GPT) and 62.9% to 34.5% (Gemini). On PBMC3k, Claude and Gemini kept exactly the
  tutorial's 2,638 cells in every successful run; on Hao no run used the tutorial thresholds.

## Run-record note
- A bug in the resume check of `task3_agent.py` (it compared run ids without the `hao_` prefix) made the first
  two Hao attempts restart from run 01, overwriting the pilot runs and the first attempt's runs. The saved runs
  come from the last attempts with unchanged agent-facing code (instruction, tools, limits); the bug is fixed.
  A network outage and a full disk also interrupted runs; both are handled now (balance checks retry; each
  run deletes the agent's intermediate .h5ad files after its results are copied).

## Post hoc condition B (SEPARATE, added after the 15-step results): 30-step cap on Hao, Claude and GPT
- Identical to the 15-step runs except the step cap (30); 5 runs per model; $1.50 per-run cap unchanged
  (`task3_agent.py --dataset hao --max-steps 30`, results in `results/task3_hao_30steps/`).

| Run | Status | Cells kept | Strict (30 types) | Lenient | Coarse (lineage) | ARI (l2) | Labels | Steps | Errors | Cost (est.) | Stop reason |
|---|---|---|---|---|---|---|---|---|---|---|---|
| hao_s30_claude_run01 | ok | 7838 | 52.8% | 67.4% | 75.6% | 0.57 | 23 | 22 | 0 | $0.705 | agent finished |
| hao_s30_gpt_run01 | ok | 7841 | 41.3% | 60.8% | 74.4% | 0.56 | 21 | 10 | 1 | $0.161 | agent finished |
| hao_s30_claude_run02 | ok | 7838 | 50.2% | 64.0% | 71.4% | 0.53 | 20 | 17 | 0 | $0.490 | agent finished |
| hao_s30_gpt_run02 | ok | 7841 | 47.5% | 65.1% | 76.5% | 0.62 | 22 | 10 | 1 | $0.146 | agent finished |
| hao_s30_claude_run03 | ok | 7841 | 49.5% | 61.0% | 66.2% | 0.53 | 21 | 26 | 1 | $1.010 | agent finished |
| hao_s30_gpt_run03 | ok | 7841 | 50.8% | 61.2% | 71.5% | 0.55 | 20 | 11 | 2 | $0.235 | agent finished |
| hao_s30_claude_run04 | ok | 7838 | 49.6% | 61.0% | 66.2% | 0.53 | 20 | 19 | 1 | $0.503 | agent finished |
| hao_s30_gpt_run04 | ok | 7841 | 35.6% | 57.4% | 72.9% | 0.51 | 18 | 8 | 0 | $0.136 | agent finished |
| hao_s30_claude_run05 | ok | 7838 | 53.2% | 64.5% | 73.1% | 0.57 | 23 | 26 | 1 | $0.963 | agent finished |
| hao_s30_gpt_run05 | ok | 7841 | 46.3% | 59.6% | 67.3% | 0.48 | 19 | 11 | 1 | $0.166 | agent finished |

| Model | Runs with labels | Strict mean, all runs (SD) | Strict, runs with labels: mean (range) | Coarse, runs with labels | ARI (l2), runs with labels | Steps mean (range) | Errors / run | Cost / run (est.) |
|---|---|---|---|---|---|---|---|---|
| Claude Sonnet 5 | 5/5 | 51.1% (1.8%) | 51.1% (49.5%-53.2%) | 70.5% | 0.55 | 22.0 (17-26) | 0.6 | $0.734 |
| GPT-5.6 Terra | 5/5 | 44.3% (6.0%) | 44.3% (35.6%-50.8%) | 72.5% | 0.54 | 10.0 (8-11) | 1.0 | $0.169 |

- **Claude:** with 30 steps all 5 runs finished (17-26 steps, mean 22), vs 3 of 7 at 15 steps, and cost rose to
  $0.734 per run. Accuracy of finished runs barely changed: 51.1% (30 steps) vs 49.7% (15 steps). Its 15-step
  failures were budget failures, not a sign that more steps give better answers.
- **GPT** used 8-11 steps with a 30-step cap, as with 15 (8-13), so the cap never bound; strict 44.3% vs 41.7%.
- **The remaining errors are decision errors:** the same T-cell confusions persist with 30 steps (CD8 T cells
  given CD4-lineage labels: Claude 55%, GPT 52%; "other T" called CD4/CD8: Claude 81%, GPT 58%).

## Post hoc condition C (SEPARATE, added after the Gemini failures): Gemini with a stateful Python tool
- Same as the main runs but `run_python` runs in one persistent sandboxed Python session per run (variables
  persist, like a notebook; same sandbox profile, 60 s per call, output truncation; a timeout restarts the
  session and says so). The tool description says variables persist. 7 runs on PBMC3k and 7 on Hao
  (`task3_agent.py [--dataset hao] --stateful --models gemini`, results in `results/task3_stateful/` and
  `results/task3_hao_stateful/`). One PBMC3k run was repeated because OpenRouter rate-limited it at step 2
  (an infrastructure failure, not an agent failure).

| Gemini 3.8 Flash | Runs with labels | Strict mean, all runs | Strict, runs with labels | Code errors (total) | Runs using all 15 steps |
|---|---|---|---|---|---|
| PBMC3k, stateless (main) | 6/7 | 53.9% | 62.9% | 26 | 7/7 |
| PBMC3k, stateful (post hoc) | 6/7 | 72.5% | 84.6% | 6 | 6/7 |
| Hao, stateless (main) | 6/7 | 29.5% | 34.5% | 18 | 7/7 |
| Hao, stateful (post hoc) | 6/7 | 21.8% | 25.5% | 0 | 7/7 |

- The stateful tool removed almost all code errors (PBMC3k 26 -> 6, Hao 18 -> 0) but Gemini still used all 15
  steps in 13 of 14 runs and still failed to save labels.csv in 1 of 7 runs per dataset.
- PBMC3k: higher accuracy (84.6% vs 62.9% for runs with labels), though one stateful run still gave coarse labels
  (42.2%). Hao: lower accuracy (25.5% vs 34.5%) with fewer distinct labels (11-18 vs 15-19) - coarser answers.
- PBMC3k stateful runs also kept exactly the tutorial's 2,638 cells.

---

# Cost (logged per answer in the results files, US$)

| Component | Cost |
|---|---|
| Task 1, PBMC3k pilot (3 models x 3 runs) | $10.74 |
| Task 1, Hao, plain models | $8.65 |
| Task 1, agents (generic, specialist incl. all 600, nudge) | $9.19 |
| Task 2, plain models | $1.80 |
| Task 2, specialist agent | $1.55 |
| Task 3, 21 agent runs (estimated at list prices) | $2.91 |
| Task 1, plain GPT and Gemini on the missing noise questions (for the 96-question comparison) | $0.24 |
| Task 1, specialist agent with GPT ($0.52) and Gemini ($0.47), 96 questions | $0.99 |
| Task 3 on Hao, 21 runs (estimated at list prices) | $4.63 |
| Post hoc A: Task 2 specialist agent with GPT ($1.09) and Gemini ($1.10) | $2.19 |
| Post hoc B: Task 3 Hao, 30-step cap, 10 runs | $4.52 |
| Post hoc C: stateful Gemini, PBMC3k ($0.38) and Hao ($0.51) | $0.89 |
| **Total logged** | **$48.31** |
| Judge check (second pass, 1,140 decisions; from the OpenRouter balance, not logged per call) | $5.63 |

Not logged: LLM-judge calls, one-question tests, and runs discarded before cost logging was added.
For the runs added in this round (the last six rows, $13.46 logged) the OpenRouter balance fell by $16.80
($26.40 -> $9.60). Most of the difference is the Task 3 Hao runs that were overwritten by the resume bug
(see "Run-record note"), plus judge calls routed through OpenRouter.

**Account totals (actual spend): [PLACEHOLDER - to be filled in from the Anthropic and OpenRouter accounts]**
