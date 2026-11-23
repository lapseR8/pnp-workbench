import numpy as np, sys, scipy.linalg as sl
from search import Prob
def iso_start(P,rng):
    D=P.D
    Q=np.zeros((D,D))
    for mon in P.mons[2]:
        a,b=mon; Q[a,b]+=1; Q[b,a]+=1
    cons=[np.ones(D)]
    V=[]
    for k in range(P.d):
        N=sl.null_space(np.array(cons))
        u=N@(rng.normal(size=N.shape[1])+1j*rng.normal(size=N.shape[1]))
        w=N@(rng.normal(size=N.shape[1])+1j*rng.normal(size=N.shape[1]))
        # Q(u+s w)=0: a s^2 + b s + c
        a=w@Q@w; b=2*u@Q@w; c=u@Q@u
        s=np.roots([a,b,c])[0]
        v=u+s*w; v/=np.linalg.norm(v)
        V.append(v); cons.append(Q@v)
    return np.array(V)
def solve(P,seed,maxn=1500,tol=1e-12):
    rng=np.random.default_rng(seed)
    A=iso_start(P,rng)
    lam=1e-2
    for it in range(maxn):
        q,_=np.linalg.qr(A.T); A=q.T
        r,J=P.F(A); nr=np.linalg.norm(r)
        if nr<tol: return A,nr,it
        JH=J.conj().T; H=JH@J; g=JH@r
        while True:
            dx=np.linalg.solve(H+lam*np.eye(H.shape[0]),-g)
            A2=A+dx.reshape(A.shape); q,_=np.linalg.qr(A2.T); A2=q.T
            r2,_=P.F(A2)
            if np.linalg.norm(r2)<nr: A=A2; lam=max(lam/3,1e-12); break
            lam*=4
            if lam>1e10: return A,nr,it
    return A,nr,maxn
if __name__=="__main__":
    m,n,d,trials=map(int,sys.argv[1:5]); s0=int(sys.argv[5]) if len(sys.argv)>5 else 0
    P=Prob(m,n,d)
    for t in range(s0,s0+trials):
        A,nr,it=solve(P,t)
        print(t,'resid',nr,'iters',it,flush=True)
        if nr<1e-9: np.save(f'sol2_{m}x{n}_d{d}_{t}.npy',A)
