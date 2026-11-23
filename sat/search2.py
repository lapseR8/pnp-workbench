import sys, itertools
from pysat.solvers import Cadical153
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool

import json
G = json.loads(sys.argv[4]) if len(sys.argv) > 4 else None
E = [tuple(e) for e in G['E']] if G else [(0,1),(1,2),(2,3),(3,4),(4,0),(1,5),(5,6),(6,7),(7,0)]
n = G['n'] if G else 8
def f(m): return any((m>>u)&1 and (m>>v)&1 for u,v in E)
# monotone circuit: enough to be true on edges, false on maximal independent sets
trues = [(1<<u)|(1<<v) for u,v in E]
falses = [m for m in range(1<<n) if not f(m) and all(f(m|(1<<i)) for i in range(n) if not (m>>i)&1)]

def solve(N, K):
    p = IDPool(); cl = []
    T = [p.id(('t',i)) for i in range(N)]  # OR?
    S = {}
    for i in range(N):
        W = n + i
        sel = []
        for a in range(W):
            for b in range(a, W):
                S[i,a,b] = p.id(('s',i,a,b)); sel.append(S[i,a,b])
        cl += CardEnc.equals(sel, 1, vpool=p, encoding=EncType.seqcounter).clauses
    cl += CardEnc.atmost(T, K, vpool=p, encoding=EncType.seqcounter).clauses
    ins = [(m, True) for m in trues] + [(m, False) for m in falses]
    for xi, (m, want) in enumerate(ins):
        def w(j):
            if j < n: return bool((m>>j)&1)
            return p.id(('v', j-n, xi))
        for i in range(N):
            v = w(n+i)
            for a in range(n+i):
                for b in range(a, n+i):
                    s = S[i,a,b]; A, B = w(a), w(b)
                    # OR case
                    def lit(z): return z
                    for t, isor in ((T[i], True), (-T[i], False)):
                        val = None
                        if isinstance(A, bool) and isinstance(B, bool):
                            val = (A or B) if isor else (A and B)
                            cl.append([-s, -t, v if val else -v])
                            continue
                        if isinstance(A, bool) or isinstance(B, bool):
                            c, x = (A, B) if isinstance(A, bool) else (B, A)
                            if isor and c or (not isor and not c):
                                cl.append([-s, -t, v if isor else -v]); continue
                            cl += [[-s, -t, -v, x], [-s, -t, v, -x]]; continue
                        if isor:
                            cl += [[-s,-t,-v,A,B],[-s,-t,v,-A],[-s,-t,v,-B]]
                        else:
                            cl += [[-s,-t,v,-A,-B],[-s,-t,-v,A],[-s,-t,-v,B]]
        out = w(n+N-1)
        cl.append([out if want else -out])
    sol = Cadical153(bootstrap_with=cl)
    if not sol.solve(): return None
    M = set(x for x in sol.get_model() if x > 0)
    gates = []
    for i in range(N):
        op = 'or' if T[i] in M else 'and'
        a, b = next((a,b) for (j,a,b),s in S.items() if j==i and s in M)
        gates.append((op, a, b))
    return gates

K = int(sys.argv[1]); lo, hi = int(sys.argv[2]), int(sys.argv[3])
for N in range(lo, hi+1):
    g = solve(N, K)
    print(N, 'SAT' if g else 'UNSAT', g, flush=True)
    if g:
        # verify on all inputs
        ok = True
        for m in range(1<<n):
            w = [bool((m>>j)&1) for j in range(n)]
            for op,a,b in g: w.append((w[a] or w[b]) if op=='or' else (w[a] and w[b]))
            ok &= (w[-1] == f(m))
        print('verified', ok); break
