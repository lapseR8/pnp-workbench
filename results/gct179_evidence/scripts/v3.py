import numpy as np, sys
from search import Prob
m,n,d,trials,s0=map(int,sys.argv[1:6])
P=Prob(m,n,d,seed=7)
D=20
# fixed: column 0 zero-sum vectors
Fx=np.zeros((3,D),complex)
for i in range(3): Fx[i,i*5+0]=1; Fx[i,3*5+0]=-1
Fx=np.linalg.qr(Fx.T)[0].T
nf=3
def F(Ap):
    A=np.vstack([Fx,Ap])
    r,J=P.F(A)
    return r,J[:,nf*D:]
for t in range(s0,s0+trials):
    rng=np.random.default_rng(t)
    Ap=rng.normal(size=(d-nf,D))+1j*rng.normal(size=(d-nf,D))
    lam=1e-2; nr=None
    for it in range(1500):
        # orthogonalize Ap against Fx and itself
        Ap=Ap-(Ap@Fx.conj().T)@Fx; Ap=np.linalg.qr(Ap.T)[0].T
        r,J=F(Ap); nr=np.linalg.norm(r)
        if nr<1e-12: break
        JH=J.conj().T; H=JH@J; g=JH@r; ok=False
        while lam<1e10:
            dx=np.linalg.solve(H+lam*np.eye(H.shape[0]),-g)
            A2=Ap+dx.reshape(Ap.shape); A2=A2-(A2@Fx.conj().T)@Fx; A2=np.linalg.qr(A2.T)[0].T
            if np.linalg.norm(F(A2)[0])<nr: Ap=A2; lam=max(lam/3,1e-12); ok=True; break
            lam*=4
        if not ok: break
    print(t,'resid',nr,'iters',it,flush=True)
    if nr<1e-9: np.save(f'v3_{m}x{n}_d{d}_{t}.npy',np.vstack([Fx,Ap]))
