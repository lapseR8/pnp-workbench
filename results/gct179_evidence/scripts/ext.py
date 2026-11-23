# exact (over Q via sympy) : vectors v with span(L6,v) valid, L6 = cross a/b
import sympy as sp, itertools
from core import *
from cross import cross
t=sp.symbols('v0:20')
def polar_sym(k,xs):
    tot=0
    for perm in itertools.permutations(range(k)):
        for mon in MON[k]:
            p=1
            for m in range(k): p*=xs[perm[m]][mon[m]]
            tot+=p
    return sp.expand(tot)
for var in ("ab" if __name__=="__main__" else ""):
    Bs=cross(variant=var); m=len(Bs)
    v=list(t)
    # linear conditions: polar_k(b_{a1..a_{k-1}}, v)=0
    lin=[]
    for k in range(1,5):
        for ms in itertools.combinations_with_replacement(range(m),k-1):
            e=polar_sym(k,[Bs[a] for a in ms]+[v])
            if e!=0: lin.append(e)
    A=sp.Matrix([[sp.diff(e,x) for x in t] for e in lin])
    ns=A.nullspace()
    print(var,'#lin eqs',len(lin),'rank',A.rank(),'solution space dim',len(ns))
    # reduce modulo L6: param v = sum c_i ns_i
    c=sp.symbols(f'c0:{len(ns)}')
    vv=sum((c[i]*ns[i] for i in range(len(ns))),sp.zeros(20,1))
    # quotient: check dim of span(ns + L6)
    M=sp.Matrix([list(n) for n in ns]+Bs)
    print('  dim span(V1+L6)=',M.rank())
    # higher conditions: polar with j copies of v, j>=2
    eqs=[]
    for k in range(2,5):
        for j in range(2,k+1):
            for ms in itertools.combinations_with_replacement(range(m),k-j):
                e=polar_sym(k,[Bs[a] for a in ms]+[list(vv)]*j)
                if e!=0: eqs.append(sp.factor(e))
    eqs=list(set(eqs))
    print('  nonlinear eqs:',len(eqs))
    for e in eqs[:30]: print('   ',e)
    G=sp.groebner(eqs,*c,order='grevlex') if eqs else None
    print('  groebner:',G)
