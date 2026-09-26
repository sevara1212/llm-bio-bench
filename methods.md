# Methods notes

Working notes for the methods section. Numbers are from `results/scores_hao.csv` unless stated.

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

- **Test set** (`src/test_judge.py`): 34 hand-written answer strings with known labels. The final judge
  scored 34/34, 33/34 and 31/34 on three repeated runs (96% overall). Residual misses are almost all
  "proliferating/cycling NK cells" mapped to NK or NK_CD56bright.
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
| Claude Sonnet 5 | 28.0% | 37.0% |
| Gemini 3.8 Flash | 22.0% | 37.7% |
| GPT-5.6 Terra | 21.7% | 32.7% |
| PanglaoDB lookup, tie-break | 16.0% | 16.0% |
| PanglaoDB lookup, fractional | 16.3% | 16.3% |
| Random guess | 3.3% | 3.3% |

(The lookups convert Ensembl IDs to symbols first, so their accuracy is the same in both formats.
Fractional numbers are over both formats together; the formats are identical for the lookups.)
Lenient: CellMarker 46.2% tie-break / 45.5% fractional; PanglaoDB 34.5% / 35.9%.
PBMC3k (8 types), strict: CellMarker 64.4% / 61.0%, PanglaoDB 36.9% / 41.9% (tie-break / fractional).

Statement of the gap: with gene symbols, the three LLMs score 4.0-9.0 points above the CellMarker
lookup (tie-break) and 5.2-10.2 points above it (fractional). With Ensembl IDs, Claude scores 0.7
points below the lookup (tie-break) and 0.5 points above it (fractional); Gemini and GPT score 5.5-7.0
points below it (5.5-5.8 fractional, 6.7-7.0 tie-break). Each LLM loses 9.0-15.7 points when the same genes are given as
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

## Open decision

- **"likely X" is currently treated as part of the primary answer** in the hand-written test cases
  (e.g. "Naive T cell (likely naive CD4+ T cell)" expected CD4 Naive), but the judge sometimes treats it
  as a hedge and answers "T cell (subtype unclear)". The hand-check rule named only "possibly" and "or".
  Decide whether "likely" is a hedge, then fix the test expectation or add it to the hedge list.
