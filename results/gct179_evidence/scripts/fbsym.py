# Symbolic V1 for F_b over P^1 (homogeneous form, polynomial in b)
import sympy as sp, itertools, random
from core import *
from ext import polar_sym
b,s=sp.symbols('b s')
def M(d):
    v=[0]*20
    for (i,j),val in d.items(): v[idx(i,j)]=val
    return v
R3=sp.Rational(1,3)
def fam(b):
    return [M({(0,0):1,(1,0):-1}),M({(3,0):1,(1,0):-1}),M({(0,2):-1,(1,2):-1,(2,2):3,(3,2):-1}),
     M({(0,1):1-b,(0,3):-b,(0,4):2*b-1}),M({(2,0):1,(2,1):-b,(2,3):b-1}),M({(2,0):-3,(2,1):1,(2,3):1,(2,4):1})]
def fam_inf():
    return [M({(0,0):1,(1,0):-1}),M({(3,0):1,(1,0):-1}),M({(0,2):-1,(1,2):-1,(2,2):3,(3,2):-1}),
     M({(0,1):-1,(0,3):-1,(0,4):2}),M({(2,1):-1,(2,3):1}),M({(2,0):-3,(2,1):1,(2,3):1,(2,4):1})]
v=list(sp.symbols('v0:20'))
def checkvalid(Bs):
    x=sp.symbols('x0:6'); t=[sp.expand(sum(x[k]*Bs[k][j] for k in range(6))) for j in range(20)]
    for k in range(1,5):
        e=sp.expand(sum(sp.prod([t[j] for j in m]) for m in MON[k]))
        if e!=0: return False
    return True
def linsys(Bs):
    lin=[]
    for k in range(1,5):
        for ms in itertools.combinations_with_replacement(range(6),k-1):
            e=polar_sym(k,[Bs[a] for a in ms]+[v])
            if e!=0: lin.append(e)
    return sp.Matrix([[sp.expand(sp.diff(e,x)) for x in v] for e in lin])
if __name__=='__main__':
    Bs=fam(b)
    print('valid symbolic in b:',checkvalid(Bs))
    print('valid at b=inf:',checkvalid(fam_inf()))
    A=linsys(Bs)
    print('A shape',A.shape,'generic rank',A.subs(b,sp.Rational(7,13)).rank())
    # rank-drop locus: gcd of determinants of random 13x13 combos
    random.seed(1)
    g=0
    for trial in range(6):
        Rm=sp.Matrix(13,A.shape[0],lambda i,j:random.randint(-5,5))
        Cm=sp.Matrix(20,13,lambda i,j:random.randint(-5,5))
        d=sp.expand((Rm*A*Cm).det(method='berkowitz'))
        g=sp.gcd(g,d)
        print('trial',trial,'gcd so far',sp.factor(g),flush=True)
    # also the 6 basis vectors must stay independent: Pluecker check
