import sympy as sp, itertools, sys
from core import *
from ext import polar_sym
def famF(b):
    a=b/(1-b); c=3*b-1
    def M(d):
        v=[0]*20
        for (i,j),val in d.items(): v[idx(i,j)]=val
        return v
    R3=sp.Rational(1,3)
    return [M({(0,0):1,(1,0):-1}),M({(3,0):1,(1,0):-1}),M({(0,2):-R3,(1,2):-R3,(2,2):1,(3,2):-R3}),
    M({(0,1):1,(0,3):-a,(0,4):a-1}),M({(2,0):1,(2,1):-b,(2,3):b-1}),M({(2,1):-c,(2,3):c-1,(2,4):1})]
def extend(Bs,label):
    m=len(Bs); t=sp.symbols('v0:20'); v=list(t)
    lin=[]
    for k in range(1,5):
        for ms in itertools.combinations_with_replacement(range(m),k-1):
            e=polar_sym(k,[Bs[a] for a in ms]+[v])
            if e!=0: lin.append(e)
    A=sp.Matrix([[sp.diff(e,x) for x in t] for e in lin])
    ns=A.nullspace(simplify=True)
    Mx=sp.Matrix([list(n) for n in ns]+Bs)
    print(label,'V1 dim',len(ns),'span(V1+L)',Mx.rank(simplify=True),flush=True)
    c=sp.symbols(f'c0:{len(ns)}')
    vv=[sp.simplify(sum(c[i]*ns[i][j] for i in range(len(ns)))) for j in range(20)]
    eqs=[]
    for k in range(2,5):
        for j in range(2,k+1):
            for ms in itertools.combinations_with_replacement(range(m),k-j):
                e=sp.factor(polar_sym(k,[Bs[a] for a in ms]+[vv]*j))
                if e!=0: eqs.append(e)
    eqs=list(set(eqs)); print('  eqs',eqs[:12])
    if eqs: print('  GB',sp.groebner(eqs,*c,order='grevlex'))
    return ns
if __name__=="__main__":
    for b in [sp.Rational(1,5),sp.Rational(2,7),sp.Rational(-3,11)]:
        Bs=famF(b)
        from modlin import valid_mod
        extend(Bs,f'F b={b}')
