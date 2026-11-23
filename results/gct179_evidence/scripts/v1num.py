import numpy as np, sys
from num import grads
from core import *
def V1(B,S=200,rng=np.random.default_rng(5)):
    m=B.shape[0]
    C=rng.normal(size=(S,m))+1j*rng.normal(size=(S,m)); X=C@B
    G=np.vstack([grads(k,X)[1] for k in range(1,5)])
    s=np.linalg.svd(G,compute_uv=False)
    _,_,vh=np.linalg.svd(G); null=vh[np.sum(s>1e-7*s[0]):].conj()
    return s, null
for f in sys.argv[1:]:
    B=np.load(f); s,null=V1(B)
    print(f,'sing vals tail',np.round(s[-8:],8),'dim V1=',null.shape[0])
    # component of V1 outside L, evaluate E2,E3,E4
    q,_=np.linalg.qr(B.T); P=q@q.conj().T
    W=null.T - P@null.T  # 20 x d
    u,sv,_=np.linalg.svd(W); w=u[:,:null.shape[0]-B.shape[0]]
    for j in range(w.shape[1]):
        x=w[:,j][None,:]
        print('  extra dir',j,[complex(np.round(E(k,x)[0],6)) for k in (1,2,3,4)])
