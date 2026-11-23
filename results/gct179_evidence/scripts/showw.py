import sympy as sp
from extF import famF, extend
from cross import cross
for Bs,lab in [(cross(variant='a'),'cross a'),(famF(sp.Rational(1,5)),'F 1/5')]:
    ns=extend(Bs,lab)
    # find vector in V1 mod L: project out L by solving
    Mx=sp.Matrix(Bs)
    for n in ns:
        aug=sp.Matrix.vstack(Mx,sp.Matrix([list(n)]))
        if aug.rank()>Mx.rank():
            print(sp.Matrix(4,5,list(n)))
