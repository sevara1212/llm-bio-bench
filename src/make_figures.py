"""Final figures -> results/figures/. Every number is computed from the score files, not typed in.

Usage: python src/make_figures.py
  (run score.py hao --skip-knob noise, score_agents.py, score_agents.py --all-nonoise and
   task2_score.py first)

  task1_accuracy_600.png        (a) Task 1 strict accuracy, 600 Hao questions (methods run on all 600)
  task1_accuracy_subset96.png   (a) companion: 96-question subset, incl. generic agent and nudge
  task1_symbol_vs_ensembl.png   (b) Task 1 symbol vs Ensembl strict accuracy per method, 600 questions
  task2_direction_matched.png   (c) Task 2 direction accuracy on the 30 expression-matched pairs
  task2_no_change_share.png     (d) Task 2 share of no_change answers on the 134 up/down questions

Colours: one fixed colour per method across all figures (reference categorical palette, fixed order);
no-LLM baselines in grays.
"""
import os
import textwrap

import matplotlib.pyplot as plt
import pandas as pd

OUT = "results/figures"
os.makedirs(OUT, exist_ok=True)
SURFACE, TEXT, TEXT_2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e0"
COLOR = {  # fixed per method in every figure
    "Claude Sonnet 5": "#2a78d6", "Gemini 3.8 Flash": "#eb6834", "GPT-5.6 Terra": "#1baf7a",
    "Specialist agent": "#eda100", "Generic agent (web)": "#e87ba4", "Specialist + nudge (post hoc)": "#008300",
}
BASELINE_GRAY, RANDOM_GRAY = "#7d7c78", "#c3c2b7"  # every no-LLM baseline / random guess, in every figure
plt.rcParams.update({"font.size": 11, "axes.edgecolor": GRID, "axes.labelcolor": TEXT_2, "xtick.color": TEXT_2,
                     "ytick.color": TEXT_2, "text.color": TEXT, "font.family": "sans-serif"})


def color_of(name, _=None):
    if name in COLOR:
        return COLOR[name]
    return RANDOM_GRAY if name.lower().startswith("random") else BASELINE_GRAY


def titles(fig, ax, title, subtitle):
    """Bold title and a wrapped gray subtitle, both anchored to the figure's left edge. Call after
    tight_layout(): it then shrinks the plot area from the top to make room for them."""
    lines = [w for line in subtitle.split("\n") for w in textwrap.wrap(line, 125)]
    h = fig.get_figheight()
    fig.text(0.012, 1 - 0.22 / h, title, ha="left", va="top", fontsize=14, fontweight="bold", color=TEXT)
    fig.text(0.012, 1 - 0.55 / h, "\n".join(lines), ha="left", va="top", fontsize=10, color=TEXT_2, linespacing=1.4)
    fig.subplots_adjust(top=1 - (0.75 + 0.2 * len(lines)) / h)


def hbar(values, title, subtitle, xlabel, path, xmax=1.0, ref=None, pct=True):
    """Horizontal bars sorted as given (top = first), value label at each bar end."""
    names = list(values)
    fig, ax = plt.subplots(figsize=(10, 0.46 * len(names) + 2.2), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    for i, n in enumerate(names):
        c = color_of(n)
        v = values[n]
        ax.barh(i, v, height=0.62, color=c, edgecolor=SURFACE, linewidth=2)
        ax.text(v + xmax * 0.01, i, f"{v:.1%}" if pct else f"{v:+.1f}", va="center", fontsize=10, color=TEXT)
    if ref:
        x, label = ref
        ax.axvline(x, color=TEXT_2, linewidth=1, alpha=0.7)
        ax.text(x + xmax * 0.008, len(names) - 0.55, label, ha="left", va="bottom", fontsize=9, color=TEXT_2)
    ax.set_yticks(range(len(names)), names)
    ax.invert_yaxis()
    ax.set_xlim(0, xmax)
    ax.set_xticks([t / 10 for t in range(0, 11) if t / 10 <= min(xmax, 1.0) + 1e-9][:: 2 if xmax > 0.6 else 1])
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0%}" if pct else f"{v:g}"))
    ax.grid(axis="x", color=GRID, linewidth=1)
    ax.set_axisbelow(True)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.tick_params(length=0)
    ax.set_xlabel(xlabel)
    fig.tight_layout()
    titles(fig, ax, title, subtitle)
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)
    print("saved", path)


# ---- Task 1, 600 Hao questions --------------------------------------------------------------------
s = pd.read_csv("results/scores_hao.csv")
a = pd.read_csv("results/agents/scores_agents_hao_all_nonoise.csv")
NAMES = {"claude-sonnet-5": "Claude Sonnet 5", "gpt-5.6-terra": "GPT-5.6 Terra", "gemini-3.8-flash": "Gemini 3.8 Flash",
         "cellmarker-lookup": "CellMarker lookup (no LLM)", "panglaodb-lookup": "PanglaoDB lookup (no LLM)"}
t1 = {NAMES[m]: d for m, d in s.groupby("model")}
t1["Specialist agent"] = a[a.condition == "specialist agent"]
assert all(len(d) == 600 for d in t1.values()), {k: len(d) for k, d in t1.items()}
# plain Claude must agree between the two score files
assert abs(a[a.condition == "plain Claude"].strict.mean() - t1["Claude Sonnet 5"].strict.mean()) < 1e-9
order = ["Specialist agent", "Claude Sonnet 5", "Gemini 3.8 Flash", "GPT-5.6 Terra",
         "CellMarker lookup (no LLM)", "PanglaoDB lookup (no LLM)"]
vals = {k: t1[k].strict.mean() for k in order}
vals["Random guess (1/30)"] = 1 / 30
hbar(vals, "Task 1: cell-type accuracy (strict)",
     "Hao 2021 PBMC, 30 cell types, same 600 questions for every method (300 symbol, 300 Ensembl), one run each.\n"
     "The generic web-search agent was run on the 96-question subset only (next figure).",
     "strict accuracy (exact cell type)", f"{OUT}/task1_accuracy_600.png", xmax=0.5)

sub = pd.read_csv("results/agents/scores_agents_hao.csv")
SUBN = {"specialist_nudge (post-hoc)": "Specialist + nudge (post hoc)", "specialist agent": "Specialist agent",
        "generic agent": "Generic agent (web)", "plain Claude": "Claude Sonnet 5", "CellMarker lookup": "CellMarker lookup (no LLM)"}
order96 = ["specialist_nudge (post-hoc)", "specialist agent", "generic agent", "plain Claude", "CellMarker lookup"]
vals96 = {SUBN[c]: sub[sub.condition == c].strict.mean() for c in order96}
vals96["Random guess (1/30)"] = 1 / 30
assert all((sub.condition == c).sum() == 96 for c in order96)
hbar(vals96, "Task 1: cell-type accuracy on the 96-question subset (strict)",
     "Fixed random subset (48 symbol / 48 Ensembl, balanced over difficulty knobs), all Claude-based except the lookup.\n"
     "The nudge was designed after the specialist's error analysis on these same questions (post hoc).",
     "strict accuracy (exact cell type)", f"{OUT}/task1_accuracy_subset96.png", xmax=0.5)

# (b) symbol vs Ensembl, dot plot with a connecting line per method
fig, ax = plt.subplots(figsize=(10, 4.9), facecolor=SURFACE)
ax.set_facecolor(SURFACE)
rows = [(k, t1[k][t1[k].format == "symbol"].strict.mean(), t1[k][t1[k].format == "ensembl"].strict.mean()) for k in order]
for i, (k, sym, ens) in enumerate(rows):
    c = color_of(k)
    ax.plot([ens, sym], [i, i], color=c, linewidth=2, solid_capstyle="round", zorder=1)
    ax.scatter([sym], [i], s=90, color=c, edgecolor=SURFACE, linewidth=2, zorder=2, marker="o")
    ax.scatter([ens], [i], s=90, color=SURFACE, edgecolor=c, linewidth=2, zorder=2, marker="o")
    gap = sym - ens
    ax.text(max(sym, ens) + 0.012, i, f"symbol {sym:.1%} | Ensembl {ens:.1%} | gap {gap * 100:+.1f} pts",
            va="center", fontsize=9.5, color=TEXT)
ax.set_yticks(range(len(rows)), [r[0] for r in rows])
ax.invert_yaxis()
ax.set_xlim(0, 0.8)
ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0%}"))
ax.grid(axis="x", color=GRID, linewidth=1)
ax.set_axisbelow(True)
for side in ("top", "right", "left"):
    ax.spines[side].set_visible(False)
ax.tick_params(length=0)
ax.set_xlabel("strict accuracy  (filled dot = gene symbols, open dot = Ensembl IDs)")
fig.tight_layout()
titles(fig, ax, "Task 1: gene symbols vs Ensembl IDs",
       "Same 600 Hao questions: each asked once with gene symbols and once with Ensembl IDs (n = 300 each). "
       "The lookups convert IDs to symbols first, so they have no gap.")
fig.savefig(f"{OUT}/task1_symbol_vs_ensembl.png", dpi=160, facecolor=SURFACE)
plt.close(fig)
print("saved", f"{OUT}/task1_symbol_vs_ensembl.png")

# ---- Task 2 --------------------------------------------------------------------------------------
t2 = pd.read_csv("results/task2_scores.csv", index_col=0)
T2N = {"specialist agent (Claude + tools)": "Specialist agent", "Claude Sonnet 5": "Claude Sonnet 5",
       "GPT-5.6 Terra": "GPT-5.6 Terra", "Gemini 3.8 Flash": "Gemini 3.8 Flash",
       "baseline: coexpr_sign": "Co-expression sign (no LLM)",
       "baseline: expression_only": "Expression level only (no LLM)",
       "baseline: coexpr_sign_|r|>0.05 (fixed in advance)": "Co-expression, |r| <= 0.05 -> no_change (no LLM)",
       "baseline: collectri": "CollecTRI (no LLM)", "baseline: always_no_change": "Always no_change (no LLM)"}
assert set(T2N) == set(t2.index), set(t2.index) ^ set(T2N)
t2 = t2.rename(index=T2N)
llms = ["Specialist agent", "Claude Sonnet 5", "Gemini 3.8 Flash", "GPT-5.6 Terra"]
bases = ["Co-expression sign (no LLM)", "Expression level only (no LLM)", "Co-expression, |r| <= 0.05 -> no_change (no LLM)",
         "CollecTRI (no LLM)", "Always no_change (no LLM)"]
hbar({k: t2.at[k, "direction_acc_matched"] for k in llms + bases},
     "Task 2: direction accuracy on expression-matched pairs",
     "Norman 2019 CRISPRa (K562). 30 up/down pairs (60 questions) matched on target expression, fixed before any model run.\n"
     "A no_change answer counts as wrong. Vertical line: a random up/down guess (50%).",
     "share of up/down questions with the correct direction", f"{OUT}/task2_direction_matched.png",
     xmax=1.0, ref=(0.5, "random up/down guess"))
hbar({k: t2.at[k, "said_no_change_on_up_down"] for k in llms + bases},
     "Task 2: how often methods answer no_change when there was an effect",
     "Share of no_change answers on the 134 questions whose true label is up or down (67 each). "
     "Lower is better; 0% means the method always commits to a direction.",
     "share answered no_change", f"{OUT}/task2_no_change_share.png", xmax=1.1)

# ---- Task 3: end-to-end analysis agents, one dot per run ------------------------------------------
if os.path.exists("results/task3/task3_scores.csv"):
    t3 = pd.read_csv("results/task3/task3_scores.csv")
    T3N = {"claude": "Claude Sonnet 5", "gpt": "GPT-5.6 Terra", "gemini": "Gemini 3.8 Flash"}
    order3 = [k for k in ["claude", "gemini", "gpt"] if k in set(t3.model_key)]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), facecolor=SURFACE, sharey=True)
    for ax, (col, xlabel) in zip(axes, [("strict", "per-cell accuracy (strict)"), ("ari_labels", "ARI vs expert clusters")]):
        ax.set_facecolor(SURFACE)
        for i, k in enumerate(order3):
            d = t3[t3.model_key == k]
            ok, bad = d[d.status == "ok"], d[d.status != "ok"]
            jitter = [(j - (len(ok) - 1) / 2) * 0.05 for j in range(len(ok))]
            ax.scatter(ok[col].fillna(0), [i + y for y in jitter], s=70, color=COLOR[T3N[k]], edgecolor=SURFACE,
                       linewidth=2, zorder=3)
            if len(bad):  # failed runs (no labels.csv) are drawn at 0 as open markers
                ax.scatter([0] * len(bad), [i] * len(bad), s=70, facecolor=SURFACE, edgecolor=COLOR[T3N[k]], linewidth=2, zorder=3)
            ax.text(1.02, i, f"median {ok[col].median():.2f}" if len(ok) else "", va="center", fontsize=9, color=TEXT_2,
                    transform=ax.get_yaxis_transform())
        ax.set_yticks(range(len(order3)), [T3N[k] for k in order3])
        if ax is axes[0]:  # the y axis is shared: invert it once, not once per panel
            ax.invert_yaxis()
        ax.set_xlim(0, 1.0)
        ax.grid(axis="x", color=GRID, linewidth=1)
        ax.set_axisbelow(True)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.tick_params(length=0)
        ax.set_xlabel(xlabel)
    n_runs = t3.groupby("model_key").size().max()
    fig.tight_layout(w_pad=6)
    titles(fig, axes[0], "Task 3: end-to-end PBMC3k analysis by an agent",
           f"Each dot is one run (up to {n_runs} per model): raw counts -> QC, normalisation, clustering, markers, labels, "
           "scored against the expert labels. Open dots at 0 = run without labels.csv; medians are over runs that "
           "produced labels.csv.")
    fig.savefig(f"{OUT}/task3_runs.png", dpi=160, facecolor=SURFACE)
    plt.close(fig)
    print("saved", f"{OUT}/task3_runs.png")
