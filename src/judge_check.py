"""Judge check: re-judge every cached Hao judge decision a second time and measure agreement.

Usage: python src/judge_check.py [--limit N] [--min-balance 4]
       python src/judge_check.py --report

The second pass uses the SAME judge (claude-sonnet-5, same instructions, label list and parsing) called
through OpenRouter, independently of the first pass. For the 1,114 decisions first made via the Anthropic
API this measures route + run-to-run agreement; for the 176 first made via OpenRouter it is a pure repeat.
New labels go to results/judge_check_hao.json only - the main cache and all scores are unchanged.
The report also rescores every Task 1 method with the second-pass labels (score change per method).
"""
import json
import os
import sys
from concurrent.futures import ThreadPoolExecutor

import pandas as pd
import requests

os.environ["JUDGE_VIA_OPENROUTER"] = "1"
import judge as J  # noqa: E402

OUT = "results/judge_check_hao.json"
cache = json.load(open(J.CACHE))
prov = json.load(open(J.PROVENANCE))
check = json.load(open(OUT)) if os.path.exists(OUT) else {}


def balance():
    for _ in range(5):
        try:
            d = requests.get("https://openrouter.ai/api/v1/credits", timeout=30,
                             headers={"Authorization": f"Bearer {os.environ['OPENROUTER_API_KEY']}"}).json()["data"]
            return d["total_credits"] - d["total_usage"]
        except Exception:  # noqa: BLE001
            __import__("time").sleep(10)
    raise RuntimeError("balance check failed")


def save():
    json.dump(check, open(OUT + ".tmp", "w"), indent=1, sort_keys=True)
    os.replace(OUT + ".tmp", OUT)


def safe_judge(a):
    try:
        return J.judge(a)
    except Exception as e:  # noqa: BLE001
        print(f"  failed {a[:50]!r}: {type(e).__name__}")
        return None


def run():
    limit = int(sys.argv[sys.argv.index("--limit") + 1]) if "--limit" in sys.argv else None
    floor = float(sys.argv[sys.argv.index("--min-balance") + 1]) if "--min-balance" in sys.argv else 4.0
    todo = sorted(a for a in cache if a not in check)[:limit]
    print(f"{len(cache)} cached decisions, {len(check)} already re-judged, {len(todo)} to do; balance ${balance():.2f}")
    for start in range(0, len(todo), 40):
        b = balance()
        if b < floor:
            print(f"STOP: balance ${b:.2f} below ${floor}")
            break
        chunk = todo[start:start + 40]
        with ThreadPoolExecutor(4) as pool:
            for a, lab in zip(chunk, pool.map(safe_judge, chunk)):
                if lab:
                    check[a] = lab
        save()
        print(f"  {len(check)}/{len(cache)} done, balance ${b:.2f}", flush=True)
    print(f"balance ${balance():.2f}")


def report():
    rows = pd.DataFrame([{"answer": a, "first": cache[a], "second": check[a], "route": prov.get(a, "anthropic")}
                         for a in cache if a in check])
    rows["agree"] = rows["first"] == rows["second"]
    lin = lambda x: J.VAGUE.get(x) or L1.get(x) or x  # noqa: E731
    rows["agree_lineage"] = rows["first"].map(lin) == rows["second"].map(lin)
    print(f"re-judged {len(rows)}/{len(cache)} decisions")
    print(rows.groupby("route")[["agree", "agree_lineage"]].agg(["mean", "size"]).round(3).to_string())
    print(f"overall: exact label agreement {rows.agree.mean():.1%}, same lineage {rows.agree_lineage.mean():.1%}")
    dis = rows[~rows.agree]
    print("\nmost common disagreements (first -> second):")
    print(dis.groupby(["first", "second"]).size().sort_values(ascending=False).head(12).to_string())
    dis.to_csv("results/judge_check_disagreements.csv", index=False)

    # rescore Task 1: answers (weighted by how often each string occurs) with first vs second labels
    from scoring import scorer
    from types import SimpleNamespace
    _, score_fn, _ = scorer("hao")
    out = []
    for p in J.results_files():
        if os.path.basename(p).startswith("scores_"):
            continue
        d = pd.read_csv(p)
        d = d[d.predicted.notna()].copy()
        if d.empty:
            continue
        d["p"] = d.predicted.astype(str)
        m1 = d.p.map(lambda a: J.rule_match(a) or cache.get(a, J.UNKNOWN))
        m2 = d.p.map(lambda a: J.rule_match(a) or check.get(a, cache.get(a, J.UNKNOWN)))
        s = lambda m: sum(score_fn(SimpleNamespace(answer=t, normalized=x)) for t, x in zip(d.answer, m)) / len(d)  # noqa: E731
        out.append({"file": os.path.basename(p), "n": len(d), "strict_first": (m1 == d.answer).mean(),
                    "strict_second": (m2 == d.answer).mean(), "lenient_first": s(m1), "lenient_second": s(m2),
                    "answers_relabelled": (m1 != m2).mean()})
    t = pd.DataFrame(out)
    t["strict_diff"] = t.strict_second - t.strict_first
    pd.set_option("display.width", 200)
    print("\nTask 1 scores with first vs second judge pass (all answers in each results file)")
    print(t.round(4).to_string(index=False))
    t.to_csv("results/judge_check_rescore.csv", index=False)

    # Task 3 on Hao: per-cell strict accuracy of every run with labels, first vs second pass
    import glob
    t3 = []
    for f in sorted(glob.glob("results/task3_hao*/runs/*/scored_cells.csv")):
        d = pd.read_csv(f, dtype=str)
        m2 = d.raw.fillna("").map(lambda a: J.rule_match(a) or check.get(a, cache.get(a, J.UNKNOWN)))
        t3.append({"run": f.split("/")[-2], "condition": f.split("/")[1], "strict_first": (d.mapped == d.l2).mean(),
                   "strict_second": (m2 == d.l2).mean()})
    t3 = pd.DataFrame(t3)
    t3["diff"] = t3.strict_second - t3.strict_first
    print("\nTask 3 on Hao, strict per run: mean and largest absolute change by condition")
    print(t3.groupby("condition").agg(runs=("run", "size"), strict_first=("strict_first", "mean"),
                                      strict_second=("strict_second", "mean"),
                                      max_abs_diff=("diff", lambda x: x.abs().max())).round(4).to_string())
    t3.to_csv("results/judge_check_task3_rescore.csv", index=False)


L1 = pd.read_csv("data/hao_labels.csv").set_index("cell_type").l1.to_dict()
if __name__ == "__main__":
    report() if "--report" in sys.argv else run()
