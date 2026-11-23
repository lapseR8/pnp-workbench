#!/usr/bin/env python3
"""watch_gaps.py STATEFILE : poll gaps.json every 60 s; exit (printing the new steps) as soon as a new
posted step with weight > 1 appears. Run in the background; relaunch after each event."""
import json, sys, time, urllib.request
URL = "https://data.pequalnp.com/data/desk/gaps.json"
state = sys.argv[1]
def fetch():
    with urllib.request.urlopen(URL, timeout=30) as r:
        d = json.load(r)
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
    new = [cur[i] for i in cur if i not in seen and (cur[i].get("weight") or 0) > 1]
    seen |= set(cur); json.dump(sorted(seen), open(state, "w"))
    if new:
        for g in new:
            print(f"NEW {g['id']} w={g.get('weight')} par={g.get('parallel')} :: {g.get('title','')[:160]}")
            print(f"   context: {str(g.get('context',''))[:300]}")
        sys.exit(0)
    time.sleep(60)
