import numpy as np, sys, itertools, glob
def rk(M,tol=1e-5):
    if M.size==0: return 0
    s=np.linalg.svd(M,compute_uv=False); return int(np.sum(s>tol*max(1,s[0])))
def inv(A,m,n):
    d=A.shape[0]
    X=A.reshape(d,m,n)
    r=[rk(X[:,i,:]) for i in range(m)]
    c=[rk(X[:,:,j]) for j in range(n)]
    p={(i,k):rk(np.concatenate([X[:,i,:],X[:,k,:]],axis=1)) for i,k in itertools.combinations(range(m),2)}
    u=[d-rk(np.delete(X,i,axis=1).reshape(d,-1)) for i in range(m)]
    v=[d-rk(np.delete(X,j,axis=2).reshape(d,-1)) for j in range(n)]
    return dict(r=r,c=c,p=p,u=u,v=v)
if __name__=="__main__":
    for f in sys.argv[1:]:
        A=np.load(f); m,n=map(int,f.split('_')[1].split('x'))
        print(f,inv(A,m,n))
