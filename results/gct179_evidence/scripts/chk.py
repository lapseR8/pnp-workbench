import sympy as sp, itertools
w=sp.Symbol('w')
om=sp.Rational(-1,2)+sp.sqrt(3)*sp.I/2
A=sp.symbols('a0:5'); B=sp.symbols('b0:5')
x=[A[l]+B[l] for l in range(5)]; y=[om*A[l]+om**2*B[l] for l in range(5)]; z=[om**2*A[l]+om*B[l] for l in range(5)]
def per34(j,s):
    cols=[l for l in range(5) if l!=j]
    tot=0
    for f in itertools.permutations(cols,3):
        tot+=(1+s*x[f[0]])*(1+s*y[f[1]])*(1+s*z[f[2]])
    return sp.expand(tot)
s=sp.Symbol('s')
def e3(v): return sum(v[i]*v[j]*v[k] for i,j,k in itertools.combinations(range(len(v)),3))
SA=sum(A);SB=sum(B)
for j in [0,3]:
    P=sp.Poly(per34(j,s),s)
    c=P.all_coeffs()[::-1]
    print(j,'c0',sp.simplify(c[0]),'c1',sp.simplify(c[1]))
    q=sum(A[l]*B[l] for l in range(5))-SA*SB+A[j]*SB+SA*B[j]-2*A[j]*B[j]
    print(' quad ratio', sp.simplify(sp.expand(c[2]/q)))
    cu=e3(A[:j]+A[j+1:])+e3(B[:j]+B[j+1:])
    print(' cubic ratio', sp.simplify(sp.expand(c[3])/cu))
