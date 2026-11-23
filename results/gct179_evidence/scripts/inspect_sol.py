import numpy as np, glob, sys
for f in sorted(glob.glob('sol_*.npy')):
    A=np.load(f)
    m,n=map(int,f.split('_')[1].split('x'))
    d=A.shape[0]
    # column-space: per-cell norm of projection
    q,_=np.linalg.qr(A.T)
    cellw=np.linalg.norm(q,axis=1).reshape(m,n)
    print(f); print(np.round(cellw,3))
    # check row r, col c support
    best=None
    for r in range(m):
        for c in range(n):
            mask=np.ones((m,n),bool); mask[r,:]=False; mask[:,c]=False
            off=np.linalg.norm(cellw[mask])
            if best is None or off<best[0]: best=(off,r,c)
    print(' best row/col support residual',best)
    # row-supported and col-supported dims
    def dimint(mask):
        # dim of L ∩ {supported on mask}: null space of coords outside mask
        out=~mask.reshape(-1)
        M=A[:,out].T
        s=np.linalg.svd(M,compute_uv=False)
        return d-np.sum(s>1e-6)
    us=[dimint(np.array([[i==r for j in range(n)] for i in range(m)])) for r in range(m)]
    vs=[dimint(np.array([[j==c for j in range(n)] for i in range(m)])) for c in range(n)]
    print(' u',us,'v',vs)
