# gen.py N K graphjson out.cnf out.map : exact-N-gate fan-in-2 monotone circuit, <=K OR gates,
# strict a<b operands, every non-output gate used, no duplicate gates.
import sys, json, pickle
from pysat.formula import IDPool, CNF
from pysat.card import CardEnc, EncType
N, K = int(sys.argv[1]), int(sys.argv[2]); G = json.loads(sys.argv[3])
n = G['n']; E = [tuple(e) for e in G['E']]
f = lambda m: any((m>>u)&1 and (m>>v)&1 for u,v in E)
ins = [((1<<u)|(1<<v), True) for u,v in E] + [(m, False) for m in range(1<<n) if not f(m) and all(f(m|(1<<i)) for i in range(n) if not (m>>i)&1)]
p = IDPool(); cnf = CNF()
T = [p.id(('t',i)) for i in range(N)]
S = {}
for i in range(N):
    sel = []
    for a in range(n+i):
        for b in range(a+1, n+i):
            S[i,a,b] = p.id(('s',i,a,b)); sel.append(S[i,a,b])
    cnf.extend(CardEnc.equals(sel, 1, vpool=p, encoding=EncType.seqcounter).clauses)
if K < N: cnf.extend(CardEnc.atmost(T, K, vpool=p, encoding=EncType.seqcounter).clauses)
for i in range(N-1):  # gate i used later
    cnf.append([s for (j,a,b),s in S.items() if j > i and (a == n+i or b == n+i)])
for (i,a,b),s in S.items():  # no duplicate gates
    for j in range(i+1, N):
        if (j,a,b) in S:
            t = S[j,a,b]; cnf.append([-s,-t,T[i],-T[j]]); cnf.append([-s,-t,-T[i],T[j]])
for xi,(m,want) in enumerate(ins):
    w = lambda j: bool((m>>j)&1) if j < n else p.id(('v', j-n, xi))
    for i in range(N):
        v = w(n+i)
        for a in range(n+i):
            for b in range(a+1, n+i):
                s = S[i,a,b]; A, B = w(a), w(b)
                for t, isor in ((T[i],True),(-T[i],False)):
                    if isinstance(A,bool) and isinstance(B,bool):
                        val = (A or B) if isor else (A and B); cnf.append([-s,-t, v if val else -v]); continue
                    if isinstance(A,bool) or isinstance(B,bool):
                        c, x = (A,B) if isinstance(A,bool) else (B,A)
                        if (isor and c) or (not isor and not c): cnf.append([-s,-t, v if isor else -v]); continue
                        cnf.extend([[-s,-t,-v,x],[-s,-t,v,-x]]); continue
                    if isor: cnf.extend([[-s,-t,-v,A,B],[-s,-t,v,-A],[-s,-t,v,-B]])
                    else: cnf.extend([[-s,-t,v,-A,-B],[-s,-t,-v,A],[-s,-t,-v,B]])
    o = w(n+N-1); cnf.append([o if want else -o])
cnf.to_file(sys.argv[4])
pickle.dump({'T':T,'S':S,'n':n,'E':E}, open(sys.argv[5],'wb'))
print(N, K, 'vars', p.top, 'clauses', len(cnf.clauses))
