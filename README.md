# llm-bio-bench

A small benchmark of how well large language models answer two everyday single-cell biology questions,
with and without tools, compared against simple non-LLM baselines. **Task 1** asks a model to name the
cell type behind a list of marker genes (30 fine-grained PBMC types from Hao et al. 2021, given as gene
symbols or Ensembl IDs, under several difficulty settings). **Task 2** asks whether switching one gene
on with CRISPR activation makes another gene go up, down or not change (Norman et al. 2019, K562
cells). Claude Sonnet 5, GPT-5.6 Terra and Gemini 3.8 Flash are compared with marker-database lookups,
co-expression and TF-network baselines, and with Claude run as a tool-using agent. All numbers below are
taken from [`methods.md`](methods.md), which documents every scoring and design choice.

## The question

Do general-purpose LLMs carry usable single-cell biology knowledge, where do they fail, and does giving
them domain tools (gene lookups, marker databases, co-expression data) help or hurt?

## Datasets

| Dataset | Use | Source |
|---|---|---|
| PBMC3k (10x Genomics, via `scanpy.datasets.pbmc3k()`) | Task 1 pilot: 8 Leiden clusters | scanpy tutorial |
| Hao et al. 2021 CITE-seq PBMC reference | Task 1: 30 expert-labelled cell types (`celltype.l2`), up to 500 cells per type | CELLxGENE, `nygc multimodal pbmc` (161,764 cells) |
| Norman et al. 2019 CRISPRa Perturb-seq | Task 2: 105 single-gene activations in K562 + 11,855 control cells | scPerturb, Zenodo DOI 10.5281/zenodo.7041849 |
| CellMarker (Hao 2021 rows excluded) and PanglaoDB | Task 1 lookup baselines and agent tools | CellMarker (downloaded by hand 2026-09-24); panglaodb.se |
| CollecTRI | Task 2 TF-target baseline and agent tool | decoupler 2.2.0 / OmniPath, snapshot 2026-09-27 |

## Tasks

- **Task 1 - cell-type annotation.** "These are the top marker genes of a cluster from human PBMCs. What
  cell type is it?" Each question comes in two gene formats (symbols / Ensembl IDs) and four difficulty
  knobs: less famous genes (marker ranks 1-10 ... 31-40), fewer genes (1, 3, 5, 10), noise (random genes
  swapped in) and a ribosomal/mitochondrial gene filter. Free-text answers are mapped to the 30 labels by
  exact matching plus a blinded LLM judge. The main comparison uses the 600 non-noise Hao questions.
- **Task 2 - perturbation direction.** "In K562 cells, gene X is activated with CRISPRa. What happens to
  the expression of gene Y? Answer up, down, or no_change." 200 questions (67 up, 67 down, 66
  no_change) built from pseudobulk differential expression, with no_change targets matched on
  expression level; 30 expression-matched up/down pairs isolate direction from expression level.

## Conditions

| Condition | Task 1 | Task 2 |
|---|---|---|
| Plain LLMs: Claude Sonnet 5 (Anthropic API), GPT-5.6 Terra and Gemini 3.8 Flash (OpenRouter) | yes | yes |
| Random guess | yes | - |
| Marker lookup, no LLM (CellMarker, PanglaoDB) | yes | - |
| Co-expression sign in control cells, expression level only, CollecTRI, always no_change (no LLM) | - | yes |
| Generic agent: Claude + web search (LangGraph ReAct, <= 5 tool calls) | 96-question subset | - |
| Specialist agent: Claude + MyGene.info + CellMarker lookup | 96 subset and all 600 | - |
| Specialist + nudge (one-line system prompt, **added post hoc** after error analysis) | 96 subset | - |
| Specialist agent: Claude + MyGene.info + control-cell co-expression + CollecTRI | - | yes |

## Headline results

### Task 1 (Hao 2021, 30 cell types, same 600 questions per method, one run each; strict = exact type)

| Method | Strict | Strict, symbols | Strict, Ensembl |
|---|---|---|---|
| Specialist agent (Claude + tools) | 41.2% | 43.3% | 39.0% |
| Claude Sonnet 5 | 32.5% | 37.0% | 28.0% |
| Gemini 3.8 Flash | 29.8% | 37.7% | 22.0% |
| GPT-5.6 Terra | 27.2% | 32.7% | 21.7% |
| CellMarker lookup (no LLM, tie-break) | 28.7% | 28.7% | 28.7% |
| PanglaoDB lookup (no LLM, tie-break) | 16.0% | 16.0% | 16.0% |
| Random guess | 3.3% | 3.3% | 3.3% |

Specialist vs plain Claude, paired: right where Claude was wrong on 59 questions, the reverse on 7
(sign test p = 2.4e-11). On the 96-question subset: generic web agent 31.2%, plain Claude 30.2%,
specialist 38.5%, specialist + nudge 42.7% (post hoc).

![Task 1 accuracy](results/figures/task1_accuracy_600.png)
![Task 1 symbol vs Ensembl](results/figures/task1_symbol_vs_ensembl.png)

### Task 2 (Norman 2019 CRISPRa, 200 questions, one run each)

| Method | Accuracy | Macro-F1 | Direction, all up/down | Direction, matched pairs | no_change on up/down |
|---|---|---|---|---|---|
| GPT-5.6 Terra | 46.0% | 0.430 | 29.9% | 31.7% | 49% |
| Gemini 3.8 Flash | 45.5% | 0.410 | 26.1% | 25.0% | 58% |
| Claude Sonnet 5 | 41.5% | 0.359 | 20.9% | 15.0% | 64% |
| Specialist agent (Claude + 3 tools) | 39.0% | 0.288 | 9.7% | 1.7% | 88% |
| Expression level only (no LLM) | 49.0% | 0.391 | 73.1% | 50.0% | 0% |
| Co-expression sign (no LLM) | 42.0% | 0.345 | 61.9% | 60.0% | 2% |
| CollecTRI (no LLM) | 33.5% | 0.176 | 0.7% | 0.0% | 99% |
| Always no_change | 33.0% | 0.165 | 0.0% | 0.0% | 100% |

![Task 2 direction on matched pairs](results/figures/task2_direction_matched.png)
![Task 2 no_change share](results/figures/task2_no_change_share.png)

## Key findings - DRAFT (to be rewritten)

1. **Domain tools help cell-type annotation.** The specialist agent reached 41.2% strict accuracy on the
   600 Hao questions vs 32.5% for plain Claude and 28.7% for the CellMarker lookup it uses as a tool.
2. **Ensembl IDs cost the plain LLMs 9.0-15.7 points; the tools remove most of that gap.** Plain models
   drop from symbols to Ensembl IDs (Gemini 37.7% -> 22.0%); the specialist drops only 4.3 points
   (43.3% -> 39.0%), and with Ensembl IDs a plain lookup (28.7%) is on par with plain Claude (28.0%).
3. **Fine subtypes stay unsolved.** CD4 T subtypes reach only 2-7% strict for plain Claude, the lookup
   and the specialist alike.
4. **On perturbation direction, no LLM beats a simple baseline.** On expression-matched up/down pairs
   the LLMs score 15.0-31.7% (they often answer no_change), while the sign of control-cell co-expression
   scores 60.0%.
5. **Tools can hurt when read the wrong way.** The Task 2 agent answered no_change on 88% of up/down
   questions and did worse than plain Claude on direction (4 vs 19 discordant questions, p = 0.003); in
   105 of its 122 wrong answers it read a weak control-cell correlation (|r| <= 0.05) as "no effect".

## Limitations

- **One run per question** for Task 1 on Hao and for Task 2 (repeat runs were only done on PBMC3k); n is
  small for some comparisons (96-question subset, 60 matched-pair questions).
- **The LLM judge is a Claude model scoring Claude and other models.** It is blind to the model and the
  true label, agreed with a 30-answer hand-check on 27/30 before its rules were tightened, scores 96% on
  a 34-case test set, and is not fully deterministic (re-judging changed 3.8% of labels, 0.3 points of
  strict accuracy).
- **Possible training-data contamination:** Hao 2021 and Norman 2019 are public and may be in the models'
  training data. Hao-derived rows were removed from CellMarker; nothing can be removed from the models.
- **Task 2 labels depend on analysis choices:** the fold change is a difference of mean log-expression
  (understates changes in lowly expressed genes), Welch's t-test with BH per activated gene, and fixed
  label thresholds. Up targets are less expressed than down targets (detection alone separates them with
  AUC 0.80), which is why direction is also reported on expression-matched pairs.
- **One model per provider, one prompt, one cell line (K562).** GPT and Gemini were called through
  OpenRouter. Agents got no instructions on how to use tools (same prompt as plain models); the nudge is
  a post hoc exception. The generic web agent was only run on the 96-question subset.
- The PBMC3k pilot is summarised in `failures.md`, not in `methods.md`.

## Future work

- **Gene pairs:** Norman 2019 has 131 two-gene activations whose genes are all also activated alone -
  predict the combination from the singles, or detect genetic interactions.
- **Task 3:** to be defined.
- Agents that are told what tool outputs mean (e.g. that weak co-expression is not evidence of no
  effect), repeat runs, and more models per provider.

## How to reproduce

Python 3.14, packages in `requirements.txt`. API keys go in `.env` (git-ignored): `ANTHROPIC_API_KEY`,
`OPENROUTER_API_KEY` (and optionally `GEMINI_API_KEY` for the direct Google API). Run from the repo root:

```bash
python3 -m venv venv && source venv/bin/activate && pip install -r requirements.txt
```

**Task 1**

```bash
python src/prepare_data.py && python src/make_questions.py pbmc      # PBMC3k pilot (downloads itself)
# Hao 2021: save https://datasets.cellxgene.cziscience.com/4078abd1-063d-4b35-aa04-324dddf8244e.h5ad
#           as data/hao2021_pbmc.h5ad (2.6 GB, git-ignored)
python src/prepare_hao.py && python src/make_questions.py hao
python src/run_eval.py anthropic/claude-sonnet-5 --dataset hao --skip-knob noise
python src/run_eval.py openai/gpt-5.6-terra --dataset hao --skip-knob noise
python src/run_eval.py google/gemini-3.8-flash --provider openrouter --dataset hao --skip-knob noise
# marker databases into data/markers/: PanglaoDB_markers_27_Mar_2020.tsv.gz and human_cell_marker.txt
python src/baseline_marker_lookup.py hao cellmarker && python src/baseline_marker_lookup.py hao panglaodb
python src/agent_subset.py
python src/run_agent.py generic && python src/run_agent.py specialist && python src/run_agent.py specialist_nudge
python src/run_agent.py specialist --questions all_nonoise
python src/judge.py && python src/score.py hao --skip-knob noise && python src/compare_baselines.py
python src/score_agents.py && python src/score_agents.py --all-nonoise
```

**Task 2**

```bash
# Norman 2019: NormanWeissman2019_filtered.h5ad from Zenodo record 7041849 into data/norman2019/ (699 MB)
python src/inspect_norman.py && python src/task2_build.py && python src/task2_baselines.py
python src/run_eval.py anthropic/claude-sonnet-5 --dataset task2
python src/run_eval.py openai/gpt-5.6-terra --dataset task2
python src/run_eval.py google/gemini-3.8-flash --provider openrouter --dataset task2
python src/run_agent_task2.py
python src/task2_score.py && python src/task2_agent_analysis.py
```

**Figures:** `python src/make_figures.py` writes `results/figures/`.

All runs are resumable (answers are saved as they come in) and `--test` flags ask one question without
saving. Seeds are fixed for every sampling step.

## Cost

API spend logged per answer in the results files (US$):

| Component | Cost |
|---|---|
| Task 1, PBMC3k pilot (3 models x 3 runs) | $10.74 |
| Task 1, Hao, plain models | $8.65 |
| Task 1, agents (generic, specialist incl. all 600, nudge) | $9.19 |
| Task 2, plain models | $1.80 |
| Task 2, specialist agent | $1.55 |
| **Total logged** | **$31.94** |

Not included: LLM-judge calls, one-question tests, and runs discarded before cost logging was added
(not measured; a few dollars).
