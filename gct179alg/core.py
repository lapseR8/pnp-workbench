import itertools, numpy as np
P,Q=4,5
N=P*Q
def idx(i,j): return i*Q+j
# monomials of E_k: list of tuples of variable indices
MON={}
for k in range(1,P+1):
    L=[]
    for S in itertools.combinations(range(P),k):
        for cols in itertools.permutations(range(Q),k):
            L.append(tuple(idx(i,c) for i,c in zip(S,cols)))
    MON[k]=L
MONA={k:np.array(v) for k,v in MON.items()}
def E(k,t):
    t=np.asarray(t)
    return np.prod(t[...,MONA[k]],axis=-1).sum(-1)
def polar(k,xs):
    """full polarization (symmetric multilinear) up to factor k!: sum over monomials and perms"""
    tot=0
    for perm in itertools.permutations(range(k)):
        prod=1
        for m in range(k):
            prod=prod*np.asarray(xs[perm[m]])[...,MONA[k][:,m]]
        tot=tot+prod.sum(-1)
    return tot
