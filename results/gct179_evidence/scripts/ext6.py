from common import *
import sympy as sp
# 6-dim example: row 0 and column 0
p=sp.symbols('p0:6')
def L6(p):
    t=[0]*20
    # row 0 off col0: entries a1..a4 sum 0
    a=[p[0],p[1],p[2],-p[0]-p[1]-p[2]]
    for j in range(1,5): t[idx[(0,j)]]=a[j-1]
    b=[p[3],p[4],p[5]]
    for i in range(1,4): t[idx[(i,0)]]=b[i-1]
    t[idx[(0,0)]]=-(p[3]+p[4]+p[5])
    return t
t=L6(p)
for k in range(1,5): print(k, Ek_sym(t,k))
v=sp.symbols('v0:20')
s=sp.Symbol('s')
tv=[t[c]+s*v[c] for c in range(20)]
lin=[]
quad=[]
for k in range(1,5):
    e=sp.Poly(Ek_sym(tv,k),s)
    cs=dict(zip([m[0] for m in e.monoms()],e.coeffs()))
    c1=sp.Poly(cs.get(1,0),*p)
    lin+=c1.coeffs()
A=sp.Matrix([[sp.diff(eq,vv) for vv in v] for eq in lin])
ns=A.nullspace()
print('linear-condition solution dim',len(ns))
for n in ns: print(list(n.T))
