"""Task 2 baselines (no LLM, no API calls) on data/task2_questions.json.

  always_no_change   predict no_change for every question
  coexpr_sign        sign of the Pearson correlation of X and Y across CONTROL cells only
                     (r > 0 -> up, r < 0 -> down); X not detected in control cells -> no_change.
                     A pure sign never says no_change, so a thresholded variant is also reported:
  coexpr_sign_|r|>0.05  |r| <= 0.05 -> no_change. The 0.05 threshold was set when this script was first
                     written, before any baseline result was seen, and has not been tuned since.
  collectri          CollecTRI TF -> target link X -> Y: weight +1 -> up, -1 -> down; no link -> no_change.
                     Coverage = share of questions where X has a CollecTRI link to Y.

Reports 3-class accuracy, macro-F1 (over up, down, no_change) and direction accuracy: on the up/down
questions only, the share predicted with the correct direction (a no_change prediction counts as wrong).
"""
import json

import pandas as pd
from sklearn.metrics import accuracy_score, confusion_matrix, f1_score

CLASSES = ["up", "down", "no_change"]
q = pd.DataFrame(json.load(open("data/task2_questions.json")))
net = pd.read_csv("data/task2/collectri_human.csv")
edge = {(s, t): w for s, t, w in zip(net.source, net.target, net.weight)}


def coexpr(r, thr=0.0):
    if r is None or pd.isna(r) or abs(r) <= thr:
        return "no_change"
    return "up" if r > 0 else "down"


q["always_no_change"] = "no_change"
q["coexpr_sign"] = [coexpr(r) for r in q.coexpr_r_ctrl]
q["coexpr_sign_|r|>0.05"] = [coexpr(r, 0.05) for r in q.coexpr_r_ctrl]
q["collectri_link"] = [(x, y) in edge for x, y in zip(q.X, q.Y)]
q["collectri"] = [("up" if edge[(x, y)] > 0 else "down") if (x, y) in edge else "no_change" for x, y in zip(q.X, q.Y)]

methods = ["always_no_change", "coexpr_sign", "coexpr_sign_|r|>0.05", "collectri"]
ud = q[q.label != "no_change"]
res = pd.DataFrame({m: {"accuracy": accuracy_score(q.label, q[m]),
                        "macro_F1": f1_score(q.label, q[m], labels=CLASSES, average="macro", zero_division=0),
                        "direction_acc_up_down": (ud[m] == ud.label).mean(),
                        "said_no_change_on_up_down": (ud[m] == "no_change").mean()}
                    for m in methods}).T
res = res.rename(index={"coexpr_sign_|r|>0.05": "coexpr_sign_|r|>0.05 (threshold fixed in advance)"})
print(f"n = {len(q)} questions: {q.label.value_counts().to_dict()}\n")
print(res.round(3).to_string())
print(f"\ncoexpr: X not detected in control cells (r undefined) on {q.coexpr_r_ctrl.isna().sum()} questions")
cov = q.collectri_link.mean()
print(f"CollecTRI coverage: {q.collectri_link.sum()}/{len(q)} questions ({cov:.1%}) have an X -> Y link; "
      f"{q[q.collectri_link].X.nunique()} distinct X; X genes with any CollecTRI targets: "
      f"{q.X.isin(set(net.source)).sum()}/{len(q)} questions")
if q.collectri_link.any():
    sub = q[q.collectri_link]
    print(f"on the covered questions: accuracy {accuracy_score(sub.label, sub.collectri):.3f}, "
          f"labels {sub.label.value_counts().to_dict()}")
for m in methods[1:]:
    print(f"\nconfusion matrix, {m} (rows = truth, columns = prediction, order {CLASSES}):")
    print(pd.DataFrame(confusion_matrix(q.label, q[m], labels=CLASSES), index=CLASSES, columns=CLASSES).to_string())
q[["id", "X", "Y", "label"] + methods + ["coexpr_r_ctrl", "collectri_link"]].to_csv("results/task2_baselines.csv", index=False)
