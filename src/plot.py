"""Step 6: accuracy vs difficulty level, one panel per knob x gene format, one line per model.

Reads results/scores.csv (run score.py first). Writes results/difficulty_curves.png.
Shaded band = 95% bootstrap CI over questions.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

KNOBS = {  # knob -> (panel title, x label, level order)
    "rank_window": ("Less famous genes", "marker ranks shown", ["1-10", "11-20", "21-30", "31-40"]),
    "n_genes": ("Fewer genes", "number of genes shown", ["1", "3", "5", "10"]),
    "noise": ("Noise injection", "random genes among 10", ["0", "2", "4", "6"]),
    "ribo_filter": ("Ribosomal/MT filter", "RPS/RPL/MT- genes removed", ["off", "on"]),
}
# Reference categorical palette, slots 1-3 in fixed order (validated all-pairs for 3 series).
COLORS = ["#2a78d6", "#eb6834", "#1baf7a"]
SURFACE, TEXT, TEXT_2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e0"

scores = pd.read_csv("results/scores.csv", dtype={"level": str})
models = sorted(scores.model.unique())
rng = np.random.default_rng(0)


def mean_ci(x, n_boot=2000):
    x = np.asarray(x)
    boots = rng.choice(x, (n_boot, len(x))).mean(axis=1)
    return x.mean(), *np.percentile(boots, [2.5, 97.5])


plt.rcParams.update({"font.size": 10, "axes.edgecolor": GRID, "axes.labelcolor": TEXT_2,
                     "xtick.color": TEXT_2, "ytick.color": TEXT_2, "text.color": TEXT})
fig, axes = plt.subplots(2, len(KNOBS), figsize=(15, 7.5), sharey=True, facecolor=SURFACE)

for col, (knob, (title, xlabel, levels)) in enumerate(KNOBS.items()):
    for row, fmt in enumerate(["symbol", "ensembl"]):
        ax = axes[row, col]
        ax.set_facecolor(SURFACE)
        x = np.arange(len(levels))
        for i, (model, color) in enumerate(zip(models, COLORS)):
            dx = (i - (len(models) - 1) / 2) * 0.06  # small sideways dodge so tied lines stay visible
            sub = scores[(scores.knob == knob) & (scores.format == fmt) & (scores.model == model)]
            stats = [mean_ci(sub[sub.level == lv].strict) for lv in levels]
            mean, lo, hi = map(np.array, zip(*stats))
            ax.fill_between(x + dx, lo, hi, color=color, alpha=0.10, linewidth=0)
            ax.plot(x + dx, mean, color=color, linewidth=2, solid_capstyle="round", solid_joinstyle="round",
                    marker="o", markersize=7, markeredgecolor=SURFACE, markeredgewidth=2, label=model)
        ax.set_xticks(x, levels)
        ax.set_ylim(-0.03, 1.05)
        ax.grid(axis="y", color=GRID, linewidth=1)
        ax.set_axisbelow(True)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.tick_params(length=0)
        if row == 0:
            ax.set_title(title, fontsize=12, fontweight="bold", loc="left", color=TEXT)
        else:
            ax.set_xlabel(xlabel)
        if col == 0:
            ax.set_ylabel(f"{'Gene symbols' if fmt == 'symbol' else 'Ensembl IDs'}\naccuracy (strict)")

handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper right", ncol=len(models), frameon=False, bbox_to_anchor=(0.99, 1.0))
fig.suptitle("Cell-type accuracy as questions get harder", x=0.01, y=0.985, ha="left", fontsize=14, fontweight="bold")
fig.text(0.01, 0.925, "PBMC3k, 8 clusters. Shaded band = 95% bootstrap CI. Noise levels pool 3 random draws.",
         color=TEXT_2, fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.90))
fig.savefig("results/difficulty_curves.png", dpi=150, facecolor=SURFACE)
print("Saved results/difficulty_curves.png")
