"""Task 2 specialist agent: paired comparison with plain Claude and classification of its wrong answers.

Usage: python src/task2_agent_analysis.py   (after the full agent run)

"Weak" co-expression = |r| <= 0.05, the same threshold as the thresholded co-expression baseline,
fixed before any model was run. r is read from the agent's own coexpression_in_control tool output
for the question's (X, Y) pair; if the agent never queried that pair, it counts as "not queried".

Wrong-answer classes:
  truth up/down, answered no_change:
    weak_coexpr_read_as_no_change     queried co-expression, |r| <= 0.05
    stronger_coexpr_still_no_change   queried co-expression, |r| > 0.05
    no_change_without_coexpr          did not query co-expression for (X, Y), or r undefined
  truth up/down, answered the opposite direction:
    wrong_dir_followed_coexpr_sign    answer matches the sign of r
    wrong_dir_against_coexpr_sign     answer opposite to the sign of r
    wrong_dir_without_coexpr
  truth no_change, answered up/down:
    false_effect_followed_coexpr_sign / false_effect_against_coexpr_sign / false_effect_without_coexpr
  no_answer                           timeout, error or unparsable reply
Writes results/agents_task2/specialist_errors_task2.csv.
"""
import json
from math import comb

import pandas as pd

q = pd.DataFrame(json.load(open("data/task2_questions.json"))).set_index("id")
agent = pd.read_csv("results/agents_task2/agent-specialist_task2.csv").set_index("id")
plain = pd.read_csv("results/claude-sonnet-5_task2.csv").set_index("id")
assert len(agent) == len(q) == len(plain), (len(agent), len(q), len(plain))
WEAK = 0.05


def sign_test(a, b):
    n, k = a + b, max(a, b)
    return min(1.0, 2 * sum(comb(n, i) for i in range(k, n + 1)) / 2 ** n) if n else 1.0


def tool_facts(qid):
    t = json.load(open(f"results/agents_task2/specialist/{qid}.json"))
    x, y = q.at[qid, "X"], q.at[qid, "Y"]
    r, used = None, []
    for m in t["messages"]:
        if m["type"] != "tool":
            continue
        used.append(m["tool_name"])
        if m["tool_name"] == "coexpression_in_control":
            try:
                out = json.loads(m["content"])
            except (json.JSONDecodeError, TypeError):
                continue
            if {out.get("gene_x"), out.get("gene_y")} == {x, y}:
                r = out.get("pearson_r_control_cells")
    return r, used


facts = {i: tool_facts(i) for i in q.index}
agent["r_seen"] = [facts[i][0] for i in agent.index]
agent["tools_used"] = [",".join(facts[i][1]) or "-" for i in agent.index]
agent["label"] = q.label

# ---- paired comparison with plain Claude ----------------------------------------------------------
a_ok = agent.predicted.reindex(q.index) == q.label
p_ok = plain.predicted.reindex(q.index) == q.label
print("PAIRED vs plain Claude, same 200 questions")
for name, ids in [("3-class, all 200", q.index), ("direction, 134 up/down", q.index[q.label != "no_change"])]:
    g = int((a_ok[ids] & ~p_ok[ids]).sum())
    l = int((~a_ok[ids] & p_ok[ids]).sum())
    print(f"   {name:<24} agent right & plain wrong: {g:>3} | plain right & agent wrong: {l:>3} | sign test p = {sign_test(g, l):.3f}")
same = (agent.predicted.reindex(q.index) == plain.predicted.reindex(q.index)).mean()
print(f"   same answer as plain Claude on {same:.0%} of questions")

# ---- how the agent used its tools -----------------------------------------------------------------
print("\nTOOL USE")
print(f"   mean tool calls {agent.n_tool_calls.mean():.2f}; queried co-expression for (X, Y) on {agent.r_seen.notna().mean():.0%}; "
      f"called CollecTRI on {agent.tools_used.str.contains('collectri').mean():.0%}; gene_info on {agent.tools_used.str.contains('gene_info').mean():.0%}")
print(f"   answer distribution: {agent.predicted.value_counts(dropna=False).to_dict()} (plain Claude: {plain.predicted.value_counts(dropna=False).to_dict()})")
seen = agent[agent.r_seen.notna()]
weak = seen.r_seen.abs() <= WEAK
print(f"   when r was weak (|r| <= {WEAK}, n = {weak.sum()}): answered no_change {(seen[weak].predicted == 'no_change').mean():.0%}")
print(f"   when r was stronger (|r| > {WEAK}, n = {(~weak).sum()}): answered no_change {(seen[~weak].predicted == 'no_change').mean():.0%}; "
      f"followed the sign of r when answering up/down {((seen[~weak].predicted == 'up') == (seen[~weak].r_seen > 0))[seen[~weak].predicted.isin(['up', 'down'])].mean():.0%}")
print(f"   share of questions with |r| <= {WEAK} among those queried: {weak.mean():.0%} "
      f"(truth on those: {seen[weak].label.value_counts().to_dict()})")


# ---- classify wrong answers -------------------------------------------------------------------------
def classify(row):
    truth, pred, r = row.label, row.predicted, row.r_seen
    if pred not in ("up", "down", "no_change") or row.finish_reason != "ok":
        return "no_answer"
    has_r = r is not None and not pd.isna(r)
    if truth in ("up", "down") and pred == "no_change":
        if not has_r:
            return "no_change_without_coexpr"
        return "weak_coexpr_read_as_no_change" if abs(r) <= WEAK else "stronger_coexpr_still_no_change"
    prefix = "wrong_dir" if truth in ("up", "down") else "false_effect"
    if not has_r or r == 0:
        return f"{prefix}_without_coexpr"
    return f"{prefix}_followed_coexpr_sign" if (pred == "up") == (r > 0) else f"{prefix}_against_coexpr_sign"


wrong = agent[agent.predicted != agent.label].copy()
wrong["class"] = wrong.apply(classify, axis=1)
wrong[["X", "Y", "label", "predicted", "r_seen", "tools_used", "class", "raw"]].to_csv(
    "results/agents_task2/specialist_errors_task2.csv")
print(f"\nWRONG ANSWERS: {len(wrong)}/200")
print(wrong["class"].value_counts().to_string())
print("\nexamples of weak_coexpr_read_as_no_change:")
ex = wrong[wrong["class"] == "weak_coexpr_read_as_no_change"].head(4)
for i, r in ex.iterrows():
    print(f"   {i}: {r.X} -> {r.Y} truth {r.label} (log2FC {q.at[i, 'log2FC']:+.2f}), r = {r.r_seen:+.3f} -> no_change")
