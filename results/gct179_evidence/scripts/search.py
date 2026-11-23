import itertools, numpy as np, sys
from math import comb
from scipy.optimize import least_squares
def setup(m,n):
    cells=[(i,j) for i in range(m) for j in range(n)]
    idx={c:k for k,c in enumerate(cells)}
    mons={}
    for k in range(1,m+1):
        L=[]
        for rows in itertools.combinations(range(m),k):
            for cols in itertools.permutations(range(n),k):
                L.append([idx[(rows[a],cols[a])] for a in range(k)])
        mons[k]=np.array(L)
    return cells,idx,mons
class Prob:
    def __init__(s,m,n,d,seed=0,ks=None):
        s.m,s.n,s.d=m,n,d
        s.cells,s.idx,s.mons=setup(m,n)
        if ks is not None: s.mons={k:s.mons[k] for k in ks}
        s.D=m*n
        rng=np.random.default_rng(seed)
        N=comb(d+m-1,m)+20
        S=rng.normal(size=(N,d))+1j*rng.normal(size=(N,d))
        s.S=S/np.linalg.norm(S,axis=1,keepdims=True)
        s.N=N
    def F(s,A):
        T=s.S@A
        out=[];J=[]
        for k,M in s.mons.items():
            P=T[:,M]
            out.append(P.prod(axis=2).sum(axis=1))
            G=np.zeros((s.N,s.D),dtype=complex)
            for a in range(k):
                oth=np.prod(np.delete(P,a,axis=2),axis=2)
                for c in range(s.D):
                    mask=(M[:,a]==c)
                    if mask.any(): G[:,c]+=oth[:,mask].sum(axis=1)
            # dres/dA[i,c] = G[s,c]*S[s,i]
            J.append((G[:,None,:]*s.S[:,:,None]).reshape(s.N,-1))
        return np.concatenate(out),np.vstack(J)
    def solve(s,seed,tol=1e-12,maxn=600):
        rng=np.random.default_rng(seed)
        A=rng.normal(size=(s.d,s.D))+1j*rng.normal(size=(s.d,s.D))
        lam=1e-3
        for it in range(maxn):
            # normalize basis (orthonormalize) to fix gauge drift
            q,_=np.linalg.qr(A.T); A=q.T
            r,J=s.F(A)
            nr=np.linalg.norm(r)
            if nr<tol: return A,nr,it
            # LM step
            JH=J.conj().T
            H=JH@J
            g=JH@r
            while True:
                dx=np.linalg.solve(H+lam*np.eye(H.shape[0]),-g)
                A2=A+dx.reshape(A.shape)
                q,_=np.linalg.qr(A2.T);A2=q.T
                r2,_=s.F(A2)
                if np.linalg.norm(r2)<nr:
                    A=A2; lam=max(lam/3,1e-12); break
                lam*=4
                if lam>1e8: return A,nr,it
        r,_=s.F(A)
        return A,np.linalg.norm(r),maxn
if __name__=="__main__":
    m,n,d,trials=map(int,sys.argv[1:5])
    P=Prob(m,n,d)
    for t in range(trials):
        A,nr,it=P.solve(t)
        print(t,'resid',nr,'iters',it,flush=True)
        if nr<1e-9:
            np.save(f'sol_{m}x{n}_d{d}_{t}.npy',A)
