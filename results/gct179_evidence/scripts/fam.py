import sympy as sp, itertools
from core import *
a,b,c=sp.symbols('a b c')
x=sp.symbols('x0:6')
def M(d):
    v=[0]*20
    for (i,j),val in d.items(): v[idx(i,j)]=val
    return v
R3=sp.Rational(1,3)
Bs=[M({(0,0):1,(1,0):-1}),M({(3,0):1,(1,0):-1}),M({(0,2):-R3,(1,2):-R3,(2,2):1,(3,2):-R3}),
    M({(0,1):1,(0,3):-a,(0,4):a-1}),M({(2,0):1,(2,1):-b,(2,3):b-1}),M({(2,1):-c,(2,3):c-1,(2,4):1})]
t=[sum(x[k]*Bs[k][v] for k in range(6)) for v in range(20)]
conds=set()
for k in range(1,5):
    e=sp.Poly(sp.expand(sum(sp.prod([t[v] for v in m]) for m in MON[k])),*x)
    for co in e.coeffs(): conds.add(sp.factor(co))
print(conds)
print(sp.solve(list(conds),[a,b,c],dict=True))
