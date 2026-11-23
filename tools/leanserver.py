#!/usr/bin/env python3
"""Warm Lean REPL server: loads each distinct import set once, then checks files in seconds.
Usage: leanserver.py  (listens on SOCK; run from anywhere)."""
import json, os, socket, subprocess, threading, time
BASE = "/Users/Ethan/Library/Application Support/Claude/scratch-workspaces/b53eb492-782f-4b89-a4cd-78a811f342fa/3a48b96d-1945-44a2-a14d-15afce52a49c/scratch-2026-11-23-4f47ca"
LEAN_DIR = BASE + "/pequalnp-proofs/lean"
REPL = BASE + "/ours/repl/.lake/build/bin/repl"
SOCK = "/tmp/claude-502/leanserver.sock"
LAKE = os.path.expanduser("~/.elan/bin/lake")

class Repl:
    def __init__(self):
        self.p = subprocess.Popen([LAKE, "env", REPL], cwd=LEAN_DIR, stdin=subprocess.PIPE,
                                  stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, bufsize=1)
        self.envs = {}
        self.lock = threading.Lock()

    def send(self, obj):
        self.p.stdin.write(json.dumps(obj) + "\n\n"); self.p.stdin.flush()
        buf = []
        while True:
            line = self.p.stdout.readline()
            if line == "":
                raise RuntimeError("repl died")
            if line.strip() == "" and buf:
                return json.loads("".join(buf))
            if line.strip():
                buf.append(line)

    def env_for(self, imports):
        key = tuple(sorted(imports))
        if key not in self.envs:
            t = time.time()
            r = self.send({"cmd": "\n".join(f"import {m}" for m in key)})
            self.envs[key] = r["env"]
            print(f"loaded {key} in {time.time()-t:.1f}s", flush=True)
        return self.envs[key]

    def check(self, text):
        lines = text.split("\n")
        imps, k = [], 0
        while k < len(lines) and (lines[k].startswith("import ") or lines[k].strip() == ""):
            if lines[k].startswith("import "):
                imps.append(lines[k].split()[1])
            k += 1
        body = "\n".join(lines[k:])
        with self.lock:
            env = self.env_for(imps)
            t = time.time()
            r = self.send({"cmd": body, "env": env})
            r["elapsed"] = round(time.time() - t, 2)
            r["line_offset"] = k
            r.pop("env", None)
            for m in r.get("messages", []):
                for key in ("pos", "endPos"):
                    if m.get(key):
                        m[key]["line"] += k
            return r

def main():
    repl = Repl()
    repl.env_for(["Mathlib", "PEqualNP.Defs"])
    if os.path.exists(SOCK):
        os.unlink(SOCK)
    s = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM); s.bind(SOCK); s.listen(8)
    print("ready on", SOCK, flush=True)
    while True:
        c, _ = s.accept()
        def handle(c=c):
            data = b""
            while not data.endswith(b"\n\n"):
                chunk = c.recv(1 << 20)
                if not chunk: break
                data += chunk
            req = json.loads(data)
            try:
                text = req.get("text") or open(req["file"]).read()
                out = repl.check(text)
            except Exception as e:
                out = {"error": str(e)}
            c.sendall(json.dumps(out).encode()); c.close()
        threading.Thread(target=handle, daemon=True).start()

if __name__ == "__main__":
    main()
