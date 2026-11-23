import sympy as sp, itertools, random, time
from sympy.polys.matrices import DomainMatrix
from core import *
from fbsym import fam, fam_inf, b
# build linear system A(b) via exact multilinear evaluation with numeric b, then interpolate (entries deg<=3 in b)
def linsys_num(Bs):
    rows=[]
    for k in range(1,5):
        for ms in itertools.combinations_with_replacement(range(6),k-1):
            vecs=[Bs[a] for a in ms]
            row=[0]*20
            for mon in MON[k]:
                # derivative wrt position mon[m] of polar: sum over choice of slot for v
                for slot in range(k):
                    others=[mon[q] for q in range(k) if q!=slot]
                    # sum over assignments of vecs to others (all perms)
                    tot=0
                    for perm in set(itertools.permutations(range(k-1))):
                        p=1
                        for q,o in enumerate(others): p*=vecs[perm[q]][o]
                        tot+=p
                    row[mon[slot]]+=tot
            rows.append(row)
    return rows
t0=time.time()
pts=list(range(2,10))
mats=[linsys_num(fam(sp.Integer(x))) for x in pts]
nr=len(mats[0])
A=sp.zeros(nr,20)
for i in range(nr):
    for j in range(20):
        A[i,j]=sp.expand(sp.interpolate([(pts[k],mats[k][i][j]) for k in range(len(pts))],b))
print('built',time.time()-t0,'max deg',max(sp.degree(e,b) if e!=0 else 0 for e in A),flush=True)
# sanity vs direct at b=7/13
Bq=fam(sp.Rational(7,13)); Aq=sp.Matrix(linsys_num(Bq))
print('interp ok:',(A.subs(b,sp.Rational(7,13))-Aq).is_zero_matrix, 'rank generic',Aq.rank(),flush=True)
dA=DomainMatrix.from_Matrix(A).convert_to(sp.ZZ[b])
random.seed(3); g=None
for trial in range(5):
    Rm=DomainMatrix.from_Matrix(sp.Matrix(13,nr,lambda i,j:random.randint(-3,3))).convert_to(sp.ZZ[b])
    Cm=DomainMatrix.from_Matrix(sp.Matrix(20,13,lambda i,j:random.randint(-3,3))).convert_to(sp.ZZ[b])
    d=(Rm*dA*Cm).det()
    d=sp.Poly(d.as_expr() if hasattr(d,'as_expr') else sp.ZZ[b].to_sympy(d),b)
    g=d if g is None else sp.gcd(g,d)
    print('trial',trial,'deg',d.degree(),'gcd',sp.factor(g.as_expr()),flush=True)
# b = infinity
Ai=sp.Matrix(linsys_num(fam_inf())); print('rank at b=inf',Ai.rank())
