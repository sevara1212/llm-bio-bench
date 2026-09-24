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

## Open decision

- **"likely X" is currently treated as part of the primary answer** in the hand-written test cases
  (e.g. "Naive T cell (likely naive CD4+ T cell)" expected CD4 Naive), but the judge sometimes treats it
  as a hedge and answers "T cell (subtype unclear)". The hand-check rule named only "possibly" and "or".
  Decide whether "likely" is a hedge, then fix the test expectation or add it to the hedge list.
