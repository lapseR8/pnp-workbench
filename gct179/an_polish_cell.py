# usage: python an_polish_cell.py <gct179 dir> <npy of bases>  (attacker B analysis, 2026-11-23)
import sys, numpy as np
sys.path.insert(0, sys.argv[1])
from core import Prob
rng=np.random.default_rng(1)
sols=list(np.load(sys.argv[2])); J=np.ones(20)
def polish(B):
    d=B.shape[0]
    for it in range(8):
        Bo=np.linalg.qr(B.T)[0].T
        Q,_=np.linalg.qr(Bo.T, mode='complete'); P=Prob(d, rng=rng, npts=400); P.U=Q.T.astype(complex)
        r,Jm=P.resid(np.zeros((d,20-d),complex),True)
        dx=np.linalg.lstsq(Jm,-r,rcond=1e-10)[0]
        B=P.basis(dx.reshape(d,20-d))
    r=P.resid(np.zeros((d,20-d),complex)); return np.linalg.qr(B.T)[0].T, np.linalg.norm(r)
for k in [2,3,6,9]:
    B,rn=polish(sols[k]); d=6
    Q,_=np.linalg.qr(B.T, mode='complete'); W=Q.T
    P=Prob(d, rng=rng, npts=400); P.U=W.astype(complex)
    r,Jm=P.resid(np.zeros((d,14),complex),True)
    Jb=Jm.reshape(Jm.shape[0], d, 14); Sm=np.vstack([Jb[:,i,:] for i in range(d)])
    s=np.linalg.svd(Sm,compute_uv=False)
    best=[]
    for rr in range(4):
        for c in range(5):
            v=J.copy(); v[rr*5+c]-=20; vp=W[d:].conj()@v
            best.append((np.linalg.norm(Sm@vp)/np.linalg.norm(vp)/s[0], rr, c))
    best.sort()
    print('sol%d polished |r|=%.1e  cell (%d,%d) rel=%.1e 2nd=%.1e | floor %.1e gap %.1e'%(k,rn,best[0][1],best[0][2],best[0][0],best[1][0],s[-1]/s[0],s[-2]/s[0]))
import itertools
INJ3=list(itertools.permutations(range(4),3))
def f34(x):  # x: 3x4 -> E1,E2,E3 (sum over injections of e_k)
    out=np.zeros(3,complex)
    for f in INJ3:
        v=[x[i,f[i]] for i in range(3)]
        out+= [sum(v), v[0]*v[1]+v[0]*v[2]+v[1]*v[2], v[0]*v[1]*v[2]]
    return out
print('---- restriction of H to the 3x4 minor avoiding the distinguished cell')
for k,(rr,cc) in zip([2,3,6,9],[(3,2),(2,0),(3,4),(1,3)]):
    B,rn=polish(sols[k])
    rows=[i for i in range(4) if i!=rr]; cols=[j for j in range(5) if j!=cc]
    mx=0; mx_other=0
    for t in range(30):
        x=((rng.normal(size=6)+1j*rng.normal(size=6))@B).reshape(4,5)
        mx=max(mx,np.abs(f34(x[np.ix_(rows,cols)])).max())
        x2=x[np.ix_([i for i in range(4) if i!=(rr+1)%4],[j for j in range(5) if j!=(cc+1)%5])]
        mx_other=max(mx_other,np.abs(f34(x2)).max())
    print('sol%d cell (%d,%d): max|E1..3 of minor| = %.1e   (control minor: %.1e)'%(k,rr,cc,mx,mx_other))
