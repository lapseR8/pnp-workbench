# usage: python an_sols6.py <gct179 dir> <npy of bases>  (attacker B analysis, 2026-11-23)
import sys, time, numpy as np
sys.path.insert(0, sys.argv[1])
from core import Prob, check_basis, forms
from iso import IsoProb
rng = np.random.default_rng(int(sys.argv[2])); n=int(sys.argv[3]); d=int(sys.argv[4]) if len(sys.argv)>4 else 6
out = sys.argv[5] if len(sys.argv)>5 else None
def tdim(B, d):
    # tangent dim of Fano scheme at span(B): chart around B itself
    Q,_ = np.linalg.qr(B.T, mode='complete'); W = Q.T  # first d rows span B
    P = Prob(d, W=W, rng=rng)
    P.U = W.astype(complex)   # chart [I X] W  -> X=0 is our space
    X0 = np.zeros((d, 20-d), complex)
    r, J = P.resid(X0, True)
    s = np.linalg.svd(J, compute_uv=False)
    return np.sum(s < 1e-7*s[0]), np.linalg.norm(r), s[-8:]
def invariants(B):
    Q,_ = np.linalg.qr(B.T); Bo = Q.T.reshape(-1,4,5)
    rr = [np.linalg.matrix_rank(Bo[:,i,:], 1e-7) for i in range(4)]
    cr = [np.linalg.matrix_rank(Bo[:,:,j], 1e-7) for j in range(5)]
    sup = (np.abs(Bo).max(axis=0) > 1e-7).sum()
    return rr, cr, sup
t0=time.time(); sols=[]
for rs in range(n):
    P = Prob(d, rng=rng); th = rng.choice([0.3,1.0,3.0])*(rng.normal(size=(d,P.m-d))+1j*rng.normal(size=(d,P.m-d)))
    th,f,it = P.solve(th, iters=200)
    if f<1e-3: th,f,it2 = P.solve(th, iters=2000, tol=1e-28)
    if f<1e-16:
        B = P.basis(th); B = np.linalg.qr(B.T)[0].T
        td, rn, s = tdim(B, d)
        print(rs, '%.1e'%f, 'tdim', td, 'inv', invariants(B), '%.0fs'%(time.time()-t0), flush=True)
        sols.append(B)
if out: np.save(out, np.array(sols))
print('nsol', len(sols), 'of', n)
