#!/usr/bin/env python3
"""watch_gaps.py STATEFILE : poll gaps.json every 60 s; exit (printing the new steps) as soon as a new
posted step with weight > 1 appears. Run in the background; relaunch after each event."""
import json, subprocess, sys, time
URL = "https://data.pequalnp.com/data/desk/gaps.json"
state = sys.argv[1]
def fetch():
    d = json.loads(subprocess.run(["curl", "-s", "--max-time", "30", URL], capture_output=True, text=True, check=True).stdout)
    g = d["gaps"]; g = g if isinstance(g, list) else list(g.values())
    return {x["id"]: x for x in g}
try:
    seen = set(json.load(open(state)))
except Exception:
    seen = set(fetch()); json.dump(sorted(seen), open(state, "w"))
while True:
    try:
        cur = fetch()
    except Exception as e:
        time.sleep(60); continue
    def good(g):
        t = str(g.get("title", "")); c = str(g.get("context", ""))
        return ((g.get("weight") or 0) >= 10 and not t.startswith("A step toward") and not t.startswith("It is false that")
                and "hypothesis of a checked reduction" not in c
                and not any(w in c.lower() for w in ("open: a hard", "hard and open", "may be false", "conjecture")))
    new = [cur[i] for i in cur if i not in seen and good(cur[i])]
    seen |= set(cur); json.dump(sorted(seen), open(state, "w"))
    if new:
        for g in new:
            print(f"NEW {g['id']} w={g.get('weight')} par={g.get('parallel')} :: {g.get('title','')[:160]}")
            print(f"   context: {str(g.get('context',''))[:300]}")
        sys.exit(0)
    time.sleep(60)
