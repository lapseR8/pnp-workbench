#!/usr/bin/env python3
"""upgen.py : witness generator for "∃ C, (FanInTwo|Monotone2) C ∧ C.size ≤ K ∧ Represents A C" rows.
SAT (CaDiCaL via python-sat) over the one-a-side rows; ops AND/OR (+ NOT when FanInTwo). Emits the Lean proof.
API: solve_stmt(statement, module, imports) -> lean text or None.   CLI: upgen.py STMT_FILE MODULE OUT"""
import re, sys, itertools
from pysat.solvers import Cadical153
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool

def parse(st):
    m = re.match(r"PEqualNP\.Represents\.RepSizeAtMost \(fun \(u : Fin (\d+)\) \(v : Fin (\d+)\) => \((\d+) : ℕ\)\.testBit \((\d+) \* u\.val \+ v\.val\)\) (\d+)$", st.strip())
    if m:
        n, mm, mask, c, K = (int(m.group(i)) for i in range(1, 6))
        return "RepSize", K, n, mm, lambda u, v: (mask >> (c * u + v)) & 1
    m = re.match(r"PEqualNP\.Represents\.RepSizeAtMost \(fun u v : Fin (\d+) => \((\d+) : ℕ\)\.testBit \((\d+) \* u\.val \+ v\.val\)\) (\d+)$", st.strip())
    if m:
        n, mask, c, K = (int(m.group(i)) for i in range(1, 5))
        return "RepSize", K, n, n, lambda u, v: (mask >> (c * u + v)) & 1
    m = re.match(r"∃ C : PEqualNP\.CircuitLowerBounds\.Circuit, PEqualNP\.MetaComplexity\.(FanInTwo|Monotone2) C ∧ C\.size ≤ (\d+) ∧ PEqualNP\.Represents\.Represents \(fun u v : Fin (\d+) => \((\d+) : ℕ\)\.testBit \((\d+) \* u\.val \+ v\.val\)\) C$", st.strip())
    if m:
        kind, K, n, mask, c = m.group(1), int(m.group(2)), int(m.group(3)), int(m.group(4)), int(m.group(5))
        return kind, K, n, n, lambda u, v: (mask >> (c * u + v)) & 1
    m = re.match(r"∃ C : PEqualNP\.CircuitLowerBounds\.Circuit, PEqualNP\.MetaComplexity\.(FanInTwo|Monotone2) C ∧ C\.size ≤ (\d+) ∧ PEqualNP\.Represents\.Represents \(fun \(u : Fin (\d+)\) \(v : Fin (\d+)\) => \((\d+) : ℕ\)\.testBit \((\d+) \* u\.val \+ v\.val\)\) C$", st.strip())
    if m:
        kind, K, n, mm, mask, c = m.group(1), int(m.group(2)), int(m.group(3)), int(m.group(4)), int(m.group(5)), int(m.group(6))
        return kind, K, n, mm, lambda u, v: (mask >> (c * u + v)) & 1
    return None

def search(kind, K, n, m, A, tmo_conf=None):
    rows = [(u, v) for u in range(n) for v in range(m)]
    I = n + m; ops = ["and", "or"] + (["not"] if kind == "FanInTwo" else [])
    for N in range(1, K + 1):
        p = IDPool(); cl = []
        X = {}; Y = {}; O = {}
        for i in range(N):
            W = I + i
            xs = [p.id(("x", i, a)) for a in range(W)]; ys = [p.id(("y", i, b)) for b in range(W)]
            os_ = [p.id(("o", i, k)) for k in range(len(ops))]
            for a in range(W): X[i, a] = xs[a]
            for b in range(W): Y[i, b] = ys[b]
            for k in range(len(ops)): O[i, k] = os_[k]
            cl += CardEnc.equals(xs, 1, vpool=p, encoding=EncType.seqcounter).clauses
            cl += CardEnc.equals(ys, 1, vpool=p, encoding=EncType.seqcounter).clauses
            cl += CardEnc.equals(os_, 1, vpool=p, encoding=EncType.seqcounter).clauses
            for a in range(W):  # symmetry: x <= y
                for b in range(a):
                    cl.append([-xs[a], -ys[b]])
        for r, (u, v) in enumerate(rows):
            def val(w):
                if w < I: return (w == u) or (w == n + v)
                return p.id(("v", w - I, r))
            for i in range(N):
                z = val(I + i); W = I + i
                for a in range(W):
                    for b in range(a, W):
                        va, vb = val(a), val(b)
                        for k, op in enumerate(ops):
                            pre = [-X[i, a], -Y[i, b], -O[i, k]]
                            if op == "not":
                                if isinstance(va, bool): cl.append(pre + [z if not va else -z])
                                else: cl += [pre + [z, va], pre + [-z, -va]]
                                continue
                            if isinstance(va, bool) and isinstance(vb, bool):
                                res = (va and vb) if op == "and" else (va or vb)
                                cl.append(pre + [z if res else -z]); continue
                            if isinstance(va, bool) or isinstance(vb, bool):
                                c_, x_ = (va, vb) if isinstance(va, bool) else (vb, va)
                                if (op == "and" and not c_) or (op == "or" and c_):
                                    cl.append(pre + [z if op == "or" else -z]); continue
                                cl += [pre + [-z, x_], pre + [z, -x_]]; continue
                            if op == "and": cl += [pre + [-z, va], pre + [-z, vb], pre + [z, -va, -vb]]
                            else: cl += [pre + [z, -va], pre + [z, -vb], pre + [-z, va, vb]]
            out = val(I + N - 1)
            cl.append([out] if A(u, v) else [-out])
        s = Cadical153(bootstrap_with=cl)
        if s.solve():
            M = set(l for l in s.get_model() if l > 0); gates = []
            for i in range(N):
                a = next(a for a in range(I + i) if X[i, a] in M); b = next(b for b in range(I + i) if Y[i, b] in M)
                op = ops[next(k for k in range(len(ops)) if O[i, k] in M)]
                gates.append((op, a, b))
            return gates, I + N - 1
    return None

def check(gates, out, n, m, A):
    for u in range(n):
        for v in range(m):
            w = [i == u or i == n + v for i in range(n + m)]
            for op, a, b in gates:
                w.append(not w[a] if op == "not" else (w[a] and w[b]) if op == "and" else (w[a] or w[b]))
            if w[out] != bool(A(u, v)): return False
    return True

def lean(st, module, imports, kind, gates, out):
    gl = ", ".join(f"⟨PEqualNP.CircuitLowerBounds.GateOp.{op}, [{a}]⟩" if op == "not" else f"⟨PEqualNP.CircuitLowerBounds.GateOp.{op}, [{a}, {b}]⟩" for op, a, b in gates)
    head = "import Mathlib\n" + "".join(f"import {x}\n" for x in imports) + "\n"
    if kind == "RepSize":
        return (head + "set_option maxRecDepth 100000 in\nset_option maxHeartbeats 0 in\n"
                f"theorem {module} : {st} := by\n  unfold PEqualNP.Represents.RepSizeAtMost\n  refine ⟨⟨[{gl}], {out}⟩, ?_, ?_, ?_⟩\n"
                "  · simp [PEqualNP.MetaComplexity.Monotone2]\n  · decide\n  · decide\n")
    return (head + "set_option maxRecDepth 100000 in\nset_option maxHeartbeats 0 in\n"
            f"theorem {module} : {st} := by\n  refine ⟨⟨[{gl}], {out}⟩, ?_, ?_, ?_⟩\n"
            f"  · simp [PEqualNP.MetaComplexity.{kind}]\n  · decide\n  · decide\n")

def solve_stmt(st, module, imports=("PEqualNP.Defs",)):
    P = parse(st)
    if not P: return None
    kind, K, n, m, A = P
    if n * m > 64: return None
    r = search(kind, K, n, m, A)
    if not r: return None
    gates, out = r
    if not check(gates, out, n, m, A): return None
    return lean(st, module, list(imports), kind, gates, out)

if __name__ == "__main__":
    st = open(sys.argv[1]).read().strip(); t = solve_stmt(st, sys.argv[2])
    if t: open(sys.argv[3], "w").write(t); print("OK")
    else: print("NO"); sys.exit(1)
