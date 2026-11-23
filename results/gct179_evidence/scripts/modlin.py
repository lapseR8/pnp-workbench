import numpy as np, itertools
from core import *
PR=1000003
def rank_mod(A,p=PR):
    A=np.array(A,dtype=np.int64)%p
    r=0; rows,cols=A.shape
    for c in range(cols):
        piv=None
        for i in range(r,rows):
            if A[i,c]: piv=i;break
        if piv is None: continue
        A[[r,piv]]=A[[piv,r]]
        inv=pow(int(A[r,c]),p-2,p)
        A[r]=(A[r]*inv)%p
        nz=np.nonzero(A[:,c])[0]
        for i in nz:
            if i!=r: A[i]=(A[i]-A[i,c]*A[r])%p
        r+=1
        if r==rows: break
    return r
def nullspace_mod(A,p=PR):
    A=np.array(A,dtype=np.int64)%p
    rows,cols=A.shape; r=0; pivs=[]
    for c in range(cols):
        piv=None
        for i in range(r,rows):
            if A[i,c]: piv=i;break
        if piv is None: continue
        A[[r,piv]]=A[[piv,r]]
        inv=pow(int(A[r,c]),p-2,p); A[r]=(A[r]*inv)%p
        for i in np.nonzero(A[:,c])[0]:
            if i!=r: A[i]=(A[i]-A[i,c]*A[r])%p
        pivs.append(c); r+=1
        if r==rows: break
    free=[c for c in range(cols) if c not in pivs]
    basis=[]
    for f in free:
        v=np.zeros(cols,dtype=np.int64); v[f]=1
        for k,c in enumerate(pivs): v[c]=(-A[k,f])%p
        basis.append(v)
    return np.array(basis)
# polarization mod p, generic: polar_mod(k, list of k integer vectors)
def polar_mod(k,xs,p=PR):
    tot=0
    for perm in itertools.permutations(range(k)):
        prod=np.ones(len(MON[k]),dtype=object)
        for m in range(k):
            prod=prod*np.array(xs[perm[m]],dtype=object)[MONA[k][:,m]]
        tot+=int(prod.sum())
    return tot%p
def valid_mod(Bs,p=PR):
    """check all polarizations vanish on span of rows of Bs"""
    m=len(Bs)
    for k in range(1,5):
        for ms in itertools.combinations_with_replacement(range(m),k):
            if polar_mod(k,[Bs[a] for a in ms],p): return False,(k,ms)
    return True,None
def tangent_system(Bs,p=PR):
    """rows: linear conditions on F (m x 20 flattened)"""
    m=len(Bs); rowsout=[]
    for k in range(1,5):
        for ms in itertools.combinations_with_replacement(range(m),k):
            row=np.zeros(m*N,dtype=object)
            for r in range(k):
                others=[Bs[a] for q,a in enumerate(ms) if q!=r]
                a=ms[r]
                # linear functional in F_a: polar(others..., e_v) for each v
                for v in range(N):
                    e=[0]*N; e[v]=1
                    row[a*N+v]+=polar_mod(k,others+[e],p)
            rowsout.append([int(x)%p for x in row])
    return np.array(rowsout,dtype=np.int64)
def tangent_dim(Bs,p=PR):
    m=len(Bs); A=tangent_system(Bs,p)
    return m*N-rank_mod(A,p)-m*m
