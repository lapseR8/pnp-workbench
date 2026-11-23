#!/usr/bin/env python3
"""leanq.py FILE.lean : check a Lean file against the warm server. Prints errors/warnings, elapsed time.
Exit 0 if no errors. Add --axioms MODULE to also print axioms."""
import json, socket, sys
SOCK = "/tmp/claude-502/leanserver.sock"
path = sys.argv[1]; text = open(path).read()
if "--axioms" in sys.argv:
    text += "\n#print axioms " + sys.argv[sys.argv.index("--axioms") + 1] + "\n"
s = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM); s.connect(SOCK)
s.sendall(json.dumps({"text": text}).encode() + b"\n\n")
data = b""
while True:
    ch = s.recv(1 << 20)
    if not ch: break
    data += ch
r = json.loads(data)
if "error" in r: print("SERVER ERROR:", r["error"]); sys.exit(2)
errs = 0
for m in r.get("messages", []):
    sev = m.get("severity"); errs += sev == "error"
    if sev == "error" or "axioms" in m.get("data", "") or "--all" in sys.argv:
        print(f"{sev} line {m['pos']['line']}: {m['data'][:600]}")
for so in r.get("sorries", []): print("sorry at line", so["pos"]["line"] + r.get("line_offset", 0))
print(f"{'OK' if errs == 0 else 'FAIL'}: {errs} errors, {r.get('elapsed')}s")
sys.exit(1 if errs else 0)
