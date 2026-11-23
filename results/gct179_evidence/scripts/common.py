import itertools, numpy as np, sympy as sp
R,C=4,5
cells=[(i,j) for i in range(R) for j in range(C)]
idx={c:n for n,c in enumerate(cells)}
# monomials of E_k: list of tuples of cell indices
mons={k:[] for k in range(1,5)}
for k in range(1,5):
    for rows in itertools.combinations(range(R),k):
        for cols in itertools.permutations(range(C),k):
            mons[k].append(tuple(idx[(rows[a],cols[a])] for a in range(k)))
def Ek(t,k):
    return sum(np.prod([t[c] for c in m]) for m in mons[k])
def Ek_sym(t,k):
    return sp.expand(sum(sp.Mul(*[t[c] for c in m]) for m in mons[k]))
def per(M):
    return sum(np.prod([M[i,f[i]] for i in range(R)]) for f in itertools.permutations(range(C),R))
