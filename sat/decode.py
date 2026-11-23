import sys, pickle
M = pickle.load(open(sys.argv[1],'rb')); lits = set()
for line in open(sys.argv[2]):
    if line.startswith('v'): lits |= {int(x) for x in line.split()[1:] if int(x) > 0}
n, E = M['n'], M['E']; N = len(M['T'])
g = [('or' if M['T'][i] in lits else 'and',) + next((a,b) for (j,a,b),s in M['S'].items() if j==i and s in lits) for i in range(N)]
f = lambda m: any((m>>u)&1 and (m>>v)&1 for u,v in E)
ok = True
for m in range(1<<n):
    w = [bool((m>>j)&1) for j in range(n)]
    for op,a,b in g: w.append((w[a] or w[b]) if op=='or' else (w[a] and w[b]))
    ok &= w[-1] == f(m)
print(g); print('ORs', sum(op=='or' for op,_,_ in g), 'verified', ok)
