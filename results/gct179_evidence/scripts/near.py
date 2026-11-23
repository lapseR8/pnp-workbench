import numpy as np, sys
from num import *
from cross import cross
import num
def gdist(A,B):
    qa,_=np.linalg.qr(A.T); qb,_=np.linalg.qr(B.T)
    s=np.linalg.svd(qa.conj().T@qb,compute_uv=False)
    return np.sqrt(np.sum(np.arccos(np.clip(s,0,1))**2))
if __name__=='__main__':
  B0=np.array(cross(variant=sys.argv[1]),dtype=complex)
  for eps in [0.05,0.2,0.5]:
    for tr in range(3):
        Bp=B0+eps*(num.rng.normal(size=B0.shape)+1j*num.rng.normal(size=B0.shape))
        B,f=solve(6,iters=400,B=Bp)
        print(eps,tr,'resid %.2e'%f,'dist to cross %.3f'%gdist(B,B0),'start dist %.3f'%gdist(Bp,B0),flush=True)
        if f<1e-9: np.save(f'near_{sys.argv[1]}_{eps}_{tr}.npy',B)
