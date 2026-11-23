import numpy as np, sys
from core import *
rng=np.random.default_rng(0)
def grads(k,X):
    # X: (S,20) -> values (S,), grads (S,20)
    M=MONA[k]; F=X[:,M]  # S x nm x k
    val=np.prod(F,axis=2).sum(1)
    G=np.zeros(X.shape,dtype=complex)
    for m in range(k):
        rest=np.prod(np.delete(F,m,axis=2),axis=2) if k>1 else np.ones(F.shape[:2])
        for col in range(1):
            pass
        np.add.at(G,(np.arange(X.shape[0])[:,None],M[:,m][None,:]),rest)
    return val,G
def residual(B,C):
    X=C@B  # S x 20
    vals=[];Js=[]
    S,m=C.shape
    for k in range(1,5):
        v,G=grads(k,X)
        scale=1.0
        vals.append(v*scale)
        # d v_s / d B[a,j] = C[s,a]*G[s,j]
        Js.append((C[:,:,None]*G[:,None,:]).reshape(S,-1)*scale)
    return np.concatenate(vals),np.vstack(Js)
def solve(m,iters=300,B=None,S=None,verbose=False):
    if S is None: S=4*m*m+40
    C=(rng.normal(size=(S,m))+1j*rng.normal(size=(S,m)))/np.sqrt(m)
    if B is None: B=rng.normal(size=(m,N))+1j*rng.normal(size=(m,N))
    lam=1e-2
    r,J=residual(B,C); f=np.linalg.norm(r)
    for it in range(iters):
        A=J.conj().T@J; g=J.conj().T@r
        step=np.linalg.solve(A+lam*np.eye(A.shape[0]),-g)
        Bn=B+step.reshape(B.shape)
        # normalize: orthonormalize rows
        q,_=np.linalg.qr(Bn.T); Bn=q.T*np.sqrt(N/ m)
        rn,Jn=residual(Bn,C); fn=np.linalg.norm(rn)
        if fn<f: B,r,J,f=Bn,rn,Jn,fn; lam=max(lam/3,1e-12)
        else: lam*=4
        if f<1e-12: break
    return B,f
if __name__=="__main__":
    m=int(sys.argv[1]); rng=np.random.default_rng(int(sys.argv[2]))
    for trial in range(int(sys.argv[3]) if len(sys.argv)>3 else 5):
        B,f=solve(m)
        print(trial,'resid',f,flush=True)
        if f<1e-9: np.save(f'sol{m}_{sys.argv[2]}_{trial}.npy',B)
