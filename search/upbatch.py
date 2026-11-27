#!/usr/bin/env python3
"""upbatch.py : run upgen over all posted ∃-Represents rows (small graphs) in parallel; write ours/prep/<ID>.lean."""
import json, os, sys, time, multiprocessing as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import upgen
PREP = os.path.join(os.path.dirname(os.path.abspath(__file__)), "prep")
def work(x):
    t = time.time()
    try:
        txt = upgen.solve_stmt(x["statement"], x["module"], tuple(x.get("imports_list") or ["PEqualNP.Defs"]))
    except Exception as e:
        return x["id"], "ERR " + repr(e)[:80], time.time() - t
    if txt:
        open(os.path.join(PREP, x["id"] + ".lean"), "w").write(txt)
        return x["id"], "WITNESS", time.time() - t
    return x["id"], "none", time.time() - t
if __name__ == "__main__":
    d = json.load(open(sys.argv[1])); g = d["gaps"]; g = g if isinstance(g, list) else list(g.values())
    todo = []
    for x in g:
        P = upgen.parse(x.get("statement", ""))
        if P and P[2] * P[3] <= 36 and not os.path.exists(os.path.join(PREP, x["id"] + ".lean")):
            todo.append(x)
    todo.sort(key=lambda x: x["id"])
    sh, ns = int(os.environ.get("SHARD", "0")), int(os.environ.get("NSHARDS", "1"))
    todo = [x for i, x in enumerate(todo) if i % ns == sh]
    todo.sort(key=lambda x: -(x.get("weight") or 0))
    print("rows", len(todo), flush=True)
    with mp.Pool(int(sys.argv[2]) if len(sys.argv) > 2 else 6) as pool:
        for i, (rid, res, secs) in enumerate(pool.imap_unordered(work, todo)):
            print(f"{i+1}/{len(todo)} {rid} {res} {secs:.0f}s", flush=True)
