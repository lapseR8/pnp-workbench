# usage: python an_obstruction.py <gct179 dir> <npy of bases>  (attacker B analysis, 2026-11-23)
import sys, numpy as np
sys.path.insert(0, sys.argv[1])
from core import Prob, forms
from ext import core
rng=np.random.default_rng(0)
sols=list(np.load(sys.argv[2]))
def wdir(B):
    d=B.shape[0]; Q,_=np.linalg.qr(B.T, mode='complete'); W=Q.T
    P=Prob(d, rng=rng); P.U=W.astype(complex)
    r,J=P.resid(np.zeros((d,20-d),complex),True)
    Jb=J.reshape(J.shape[0], d, 20-d); S=np.vstack([Jb[:,i,:] for i in range(d)])
    u,s,vh=np.linalg.svd(S); return vh[-1].conj()@W[d:], s[-2]/s[0], s[-1]/s[0]
for name,B in [('rowcol',core('rowcol')),('g180',core('g180'))]+[('sol%d'%i,b) for i,b in enumerate(sols)]:
    w,s2,s1=wdir(B)
    Bo=np.linalg.qr(B.T)[0].T
    h=(rng.normal(size=(40,Bo.shape[0]))+1j*rng.normal(size=(40,Bo.shape[0])))@Bo
    ss=np.exp(2j*np.pi*np.arange(5)/5)
    vals=np.array([forms((h+s*w).reshape(-1,4,5)) for s in ss])
    coef=np.fft.fft(vals,axis=0)/5          # coef[j][n,k]: s^j coefficient of E_k(h_n + s w)
    Ew=forms(w.reshape(1,4,5))[0]
    # mixed: coefficient of s^j in E_k for k>j (should equal 0) ; s^k in E_k equals E_k(w)
    mixed=max(np.abs(coef[j][:,k]).max() for k in range(4) for j in range(1,5) if j!=k+1)
    pure=max(np.abs(coef[k+1][:,k]-Ew[k]).max() for k in range(4))
    E1,E2,E3,E4=Ew
    print('%-7s sv gap %.1e/%.1e  mixed=%.1e pure-dev=%.1e  E1(w)=%.1e  E3^2/(E2E4)=%s  E2^2/E4=%s'%(name,s2,s1,mixed,pure,abs(E1),np.round(E3**2/(E2*E4),6),np.round(E2**2/E4,5)))
print('---- is w = J - 20 e_rc mod H ?')
J=np.ones(20)
for name,B in [('rowcol',core('rowcol')),('g180',core('g180'))]+[('sol%d'%i,b) for i,b in enumerate(sols)]:
    w,_,_=wdir(B); Bo=np.linalg.qr(B.T)[0].T
    res=[]
    for r in range(4):
        for c in range(5):
            v=J.copy(); v[r*5+c]-=20
            M=np.vstack([Bo,w/np.linalg.norm(w),v/np.linalg.norm(v)])
            s=np.linalg.svd(M,compute_uv=False); res.append((s[-1],r,c))
    res.sort()
    # also: does H contain J? does H+w contain J?
    sJ=np.linalg.svd(np.vstack([Bo,w/np.linalg.norm(w),J/np.sqrt(20)]),compute_uv=False)[-1]
    print(name,'best (r,c)=(%d,%d) sv=%.1e next=%.1e | J in H+w? sv=%.1e'%(res[0][1],res[0][2],res[0][0],res[1][0],sJ))
print('---- first-order residual of v=J-20e_rc (relative), vs sv floor/gap')
for name,B in [('rowcol',core('rowcol')),('g180',core('g180'))]+[('sol%d'%i,b) for i,b in enumerate(sols)]:
    d=B.shape[0]; Q,_=np.linalg.qr(B.T, mode='complete'); W=Q.T
    P=Prob(d, rng=rng); P.U=W.astype(complex)
    r,Jm=P.resid(np.zeros((d,20-d),complex),True)
    Jb=Jm.reshape(Jm.shape[0], d, 20-d); Sm=np.vstack([Jb[:,i,:] for i in range(d)])
    s=np.linalg.svd(Sm,compute_uv=False)
    best=[]
    for rr in range(4):
        for c in range(5):
            v=J.copy(); v[rr*5+c]-=20; vp=W[d:].conj()@v
            best.append((np.linalg.norm(Sm@vp)/np.linalg.norm(vp)/s[0], rr, c))
    best.sort()
    print(name, 'cell (%d,%d) rel=%.1e  2nd=%.1e | sv floor %.1e gap %.1e'%(best[0][1],best[0][2],best[0][0],best[1][0],s[-1]/s[0],s[-2]/s[0]))
