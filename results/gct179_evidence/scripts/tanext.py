import numpy as np
from search import Prob
P=Prob(4,5,1)
def grads(t):
    out=[]
    for k,M in P.mons.items():
        G=np.zeros(20,complex)
        Pv=t[M]
        for a in range(k):
            oth=np.prod(np.delete(Pv,a,axis=1),axis=1)
            np.add.at(G,M[:,a],oth)
        out.append(G)
    return out
def T(A,ns=200,seed=0):
    rng=np.random.default_rng(seed); rows=[]
    for _ in range(ns):
        c=rng.normal(size=A.shape[0])+1j*rng.normal(size=A.shape[0])
        rows+=grads(c@A)
    M=np.array(rows); s=np.linalg.svd(M,compute_uv=False)
    return 20-np.sum(s>1e-8*s[0]), M
if __name__=="__main__":
    z=np.exp(2j*np.pi*np.arange(5)/5)
    A=np.zeros((4,20),complex)
    for i in range(4): A[i,i*5:(i+1)*5]=z
    dT,M=T(A); print('dim T(C^4 x zeta) =',dT)
    U,s,Vh=np.linalg.svd(M); N=Vh[np.sum(s>1e-8*s[0]):].conj()
    # express N in Fourier basis per row
    F=np.array([[np.exp(-2j*np.pi*j*k/5) for j in range(5)] for k in range(5)])/5
    for v in N:
        X=v.reshape(4,5)@F.T
        pass
    # which Fourier modes appear in T: project
    Xs=np.array([ (v.reshape(4,5)@F.T) for v in N])  # coefficient of chi_k in each row
    print('mode weights', np.round(np.linalg.norm(Xs,axis=(0,1)),3))
