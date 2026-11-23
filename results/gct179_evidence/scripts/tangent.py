import numpy as np, glob
from search import Prob
for f in sorted(glob.glob('sol_*.npy')):
    A=np.load(f); m,n=map(int,f.split('_')[1].split('x')); d=A.shape[0]
    P=Prob(m,n,d,seed=5)
    r,J=P.F(A)
    s=np.linalg.svd(J,compute_uv=False)
    null=J.shape[1]-np.sum(s>1e-7*s[0])
    print(f,'resid',np.linalg.norm(r),'jac nullity',null,'gauge',d*d,'=> local family dim (Grassmannian)',null-d*d, ' small svals', s[-null-3:-null+1] if null>0 else s[-3:])
