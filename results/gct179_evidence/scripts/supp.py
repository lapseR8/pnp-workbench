import numpy as np, sys
import num
from num import residual
from core import *
def solve_mask(m,mask,iters=400,B=None,rng=None):
    S=4*m*m+40
    C=(rng.normal(size=(S,m))+1j*rng.normal(size=(S,m)))/np.sqrt(m)
    if B is None: B=(rng.normal(size=(m,N))+1j*rng.normal(size=(m,N)))*mask
    sel=np.tile(mask,m).astype(bool)
    lam=1e-2; r,J=residual(B,C); f=np.linalg.norm(r)
    for it in range(iters):
        Js=J[:,sel]
        A=Js.conj().T@Js; g=Js.conj().T@r
        st=np.zeros(m*N,dtype=complex); st[sel]=np.linalg.solve(A+lam*np.eye(A.shape[0]),-g)
        Bn=B+st.reshape(B.shape)
        q,_=np.linalg.qr(Bn.T); Bn=(q.T*np.sqrt(N/m))*mask
        rn,Jn=residual(Bn,C); fn=np.linalg.norm(rn)
        if fn<f: B,r,J,f=Bn,rn,Jn,fn; lam=max(lam/3,1e-12)
        else: lam*=4
        if f<1e-12: break
    return B,f
if __name__=="__main__":
    m=int(sys.argv[1]); rows=[int(x) for x in sys.argv[2].split(',') if x]; cols=[int(x) for x in sys.argv[3].split(',') if x]
    mask=np.array([1.0 if (k//5 in rows or k%5 in cols) else 0.0 for k in range(20)])
    rng=np.random.default_rng(int(sys.argv[4]))
    for tr in range(int(sys.argv[5])):
        B,f=solve_mask(m,mask,rng=rng)
        print(m,rows,cols,tr,'resid %.2e'%f,flush=True)
        if f<1e-9: np.save(f'supp_{m}_{sys.argv[2]}_{sys.argv[3]}_{sys.argv[4]}_{tr}.npy',B)
