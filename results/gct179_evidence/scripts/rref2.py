import numpy as np, sys
from scipy.linalg import qr
np.set_printoptions(linewidth=220,precision=4,suppress=True)
B=np.load(sys.argv[1])
_,_,perm=qr(B,pivoting=True)
piv=sorted(perm[:B.shape[0]])
R=np.linalg.solve(B[:,piv],B)
print('pivots',piv)
for row in R: print(np.round(row.reshape(4,5),4)); print()
