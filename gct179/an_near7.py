# usage: python an_near7.py <gct179 dir> <npy of bases>  (attacker B analysis, 2026-11-23)
import sys, numpy as np
sys.path.insert(0, sys.argv[1])
from core import Prob, check_basis, forms
from ext import core
rng=np.random.default_rng(3)
Bs=np.load(sys.argv[2])
def split(B, npts=400):
    Bo=np.linalg.qr(B.T)[0].T; d=Bo.shape[0]
    c=rng.normal(size=(npts,d))+1j*rng.normal(size=(npts,d)); c/=np.linalg.norm(c,axis=1,keepdims=True)
    E=forms((c@Bo).reshape(-1,4,5)); return np.sqrt((np.abs(E)**2).mean(axis=0))
def best_hyper(B, n=6):
    Bo=np.linalg.qr(B.T)[0].T
    best=(1e9,None)
    for r in range(n):
        P=Prob(6, W=Bo, rng=rng); th=0.5*(rng.normal(size=(6,1))+1j*rng.normal(size=(6,1)))
        th,f,it=P.solve(th,iters=300)
        if f<best[0]: best=(f,P.basis(th))
    return best
def nearest6(H, n=1):
    # LM in full space started at H (chart centered at H)
    Q,_=np.linalg.qr(H.T, mode='complete'); P=Prob(6, rng=rng); P.U=Q.T.astype(complex)
    X,f,it=P.solve(np.zeros((6,14),complex),iters=2000,tol=1e-28)
    H6=P.basis(X); return f, H6, np.linalg.norm(X)
def principal_angles(A,B):
    qa=np.linalg.qr(A.T)[0]; qb=np.linalg.qr(B.T)[0]
    s=np.linalg.svd(qa.conj().T@qb,compute_uv=False); return np.sqrt(np.maximum(0,1-np.minimum(s,1)**2))
fams={'rowcol':core('rowcol'),'g180':core('g180')}
for k,B in enumerate(Bs[:5]):
    sp=split(B)
    fh,H=best_hyper(B)
    f6,H6,dx=nearest6(np.linalg.qr(H.T)[0].T)
    Bo=np.linalg.qr(B.T)[0].T
    # how far is the exact 6-sol from the 7-space: sin of principal angles
    ang=principal_angles(H6,Bo)
    sp6=split(H)
    print('min%d  rms E1..E4 on 7-space: %s | best hyperplane f=%.2e rms=%s | nearest exact 6-sol f=%.1e |X|=%.2e maxangle(6sol,7space)=%.2e' % (
        k, ' '.join('%.1e'%x for x in sp), fh, ' '.join('%.1e'%x for x in sp6), f6, dx, ang.max()), flush=True)
    rr=[np.linalg.matrix_rank(np.linalg.qr(H6.T)[0].T.reshape(-1,4,5)[:,i,:],1e-6) for i in range(4)]
    cr=[np.linalg.matrix_rank(np.linalg.qr(H6.T)[0].T.reshape(-1,4,5)[:,:,j],1e-6) for j in range(5)]
    sup=(np.abs(np.linalg.qr(H6.T)[0].T.reshape(-1,4,5)).max(axis=0)>1e-6).sum()
    print('   6-sol row ranks',rr,'col ranks',cr,'support',sup)
print('---- re-chart polish (orthonormal centered chart, normalized forms)')
for k,B in enumerate(Bs[:5]):
    Bo=np.linalg.qr(B.T)[0].T
    Q,_=np.linalg.qr(Bo.T, mode='complete'); P=Prob(7, rng=rng); P.U=Q.T.astype(complex)
    X0=np.zeros((7,13),complex); r0=P.resid(X0); f0=np.vdot(r0,r0).real
    X,f,it=P.solve(X0,iters=4000,tol=1e-30)
    r,J=P.resid(X,True); s=np.linalg.svd(J,compute_uv=False)
    g=np.linalg.norm(J.conj().T@r)
    print('min%d f(orthonormal start)=%.3e -> f=%.3e after %d it, |grad|=%.1e, |X|=%.2e, smallest sv of J: %s'%(k,f0,f,it,g,np.linalg.norm(X),' '.join('%.1e'%x for x in s[-4:])), flush=True)
