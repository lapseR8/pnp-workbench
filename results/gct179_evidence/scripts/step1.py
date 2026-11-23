import sympy as sp, itertools
from core import *
t=sp.symbols('t0:20')
def Ek(k): return sum(sp.prod([t[v] for v in m]) for m in MON[k])
E2=sp.expand(Ek(2))
# quadratic form matrix
B=sp.hessian(E2,t)
# restrict to H: basis of H
Hb=[]
for v in range(1,20):
    e=[0]*20; e[0]=-1; e[v]=1; Hb.append(e)
Hm=sp.Matrix(Hb).T
BH=Hm.T*B*Hm
print('rank B on C^20',B.rank(),' rank on H',BH.rank())
print('eigs B', {k:v for k,v in B.eigenvals().items()})
