import numpy as np, sys, scipy.linalg as sl
from scipy.optimize import least_squares
from search import Prob
m,n,d,trials,s0=map(int,sys.argv[1:6])
P=Prob(m,n,d,seed=3)
def build(z):
    z=z[:len(z)//2]+1j*z[len(z)//2:]
    phi=z[:5]; Ap=z[5:].reshape(d-3,20)
    # U = 1^perp ∩ phi^perp in C^5 (bilinear)
    N=sl.null_space(np.vstack([np.ones(5),phi]))  # 5x3 (uses conj? null_space solves M x=0 exactly)
    Ub=np.zeros((3,20),complex); Ub[:,0:5]=N.T
    A=np.vstack([Ub,Ap])
    return A
def res(z):
    A=build(z)
    # normalize rows of Ap block to avoid collapse: orthonormalize whole basis
    q,_=np.linalg.qr(A.T); A=q.T
    r,_=P.F(A)
    # penalize degeneracy via orthonormalization only
    return np.concatenate([r.real,r.imag])
for t in range(s0,s0+trials):
    rng=np.random.default_rng(t)
    z0=rng.normal(size=2*(5+(d-3)*20))
    out=least_squares(res,z0,method='lm',max_nfev=20000,xtol=1e-15,ftol=1e-15,gtol=1e-15)
    nr=np.linalg.norm(out.fun)
    print(t,'resid',nr,'nfev',out.nfev,flush=True)
    if nr<1e-8:
        A=build(out.x); np.save(f'u3_{m}x{n}_d{d}_{t}.npy',A)
