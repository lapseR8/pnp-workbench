import numpy as np, sys
np.set_printoptions(linewidth=200,precision=4,suppress=True)
B=np.load(sys.argv[1])
def rref(A,tol=1e-8):
    A=A.copy(); r=0; piv=[]
    for c in range(A.shape[1]):
        i=np.argmax(np.abs(A[r:,c]))+r
        if abs(A[i,c])<tol: continue
        A[[r,i]]=A[[i,r]]; A[r]/=A[r,c]
        for k in range(A.shape[0]):
            if k!=r: A[k]-=A[k,c]*A[r]
        piv.append(c); r+=1
        if r==A.shape[0]: break
    return A,piv
order=[int(x) for x in sys.argv[2].split(',')] if len(sys.argv)>2 else list(range(20))
R,piv=rref(B[:,order])
print('pivots',[order[p] for p in piv])
for row in R:
    out=np.zeros(20,dtype=complex); out[order]=row
    print(np.round(out.reshape(4,5),4))
