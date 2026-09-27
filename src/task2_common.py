"""Shared by the Task 2 plain-model runs (run_eval.py --dataset task2) and the Task 2 agent."""
import json
import re

CLASSES = ("up", "down", "no_change")


def parse_direction(text):
    """(direction, confidence) from a JSON reply; direction is the first up/down/no_change value found."""
    try:
        obj = json.loads(re.search(r"\{.*\}", text or "", re.S).group())
    except (AttributeError, json.JSONDecodeError):
        return None, None
    if not isinstance(obj, dict):
        return None, None
    conf = next((v for k, v in obj.items() if "conf" in k.lower()), None)
    for v in obj.values():
        v = str(v).strip().lower().replace(" ", "_").replace("-", "_")
        if v in CLASSES:
            return v, conf
    return None, conf
