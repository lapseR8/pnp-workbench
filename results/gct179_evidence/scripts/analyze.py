import numpy as np, sys
from near import gdist
from cross import cross
from core import *
B=np.load(sys.argv[1])
best=min(((gdist(B,np.array(cross(r,c,v),dtype=complex)),r,c,v) for r in range(4) for c in range(5) for v in 'ab'))
print('closest cross',best)
# projector onto L, diag => which coordinates involved
q,_=np.linalg.qr(B.T); Pm=q@q.conj().T
print('diag of projector (4x5):'); print(np.round(np.real(np.diag(Pm)).reshape(4,5),3))
# u_i: dim of L intersect row-i-supported
for i in range(4):
    mask=np.array([ (k//5)!=i for k in range(20)])
    # vectors in L with zero entries off row i: B^T c restricted to mask=0
    s=np.linalg.svd(B[:,mask].T,compute_uv=False)
    print('row',i,'u=',int(np.sum(s<1e-8)), ' r=',np.linalg.matrix_rank(B[:,~mask],1e-8))
for j in range(5):
    mask=np.array([ (k%5)!=j for k in range(20)])
    s=np.linalg.svd(B[:,mask].T,compute_uv=False)
    print('col',j,'u=',int(np.sum(s<1e-8)),' r=',np.linalg.matrix_rank(B[:,~mask],1e-8))
