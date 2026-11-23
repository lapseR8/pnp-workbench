import sympy as sp, itertools
def e3(v):
    return sum(v[i]*v[j]*v[k] for i,j,k in itertools.combinations(range(len(v)),3))
c,p,q,r=sp.symbols('c p q r')
cases={}
# s=2
a=[c,c,p,q,-p-q]; b=[r,-r,0,0,0]; cases[2]=(a,b,[c,p,q,r])
a=[c,c,c,p,-c-p]; b=[q,r,-q-r,0,0]; cases[3]=(a,b,[c,p,q,r])
a=[c,c,c,c,-2*c]; b=[p,q,r,-p-q-r,0]; cases[4]=(a,b,[c,p,q,r])
n1,n2,n3=sp.symbols('n1 n2 n3')
for s,(a,b,vars_) in cases.items():
    F=[sp.expand(e3(a[:j]+a[j+1:])+e3(b[:j]+b[j+1:])) for j in range(5)]
    sols=[]
    for k in range(4):
        others=[v for v in vars_ if v!=vars_[k]]
        sub={vars_[k]: n1*others[0]+n2*others[1]+n3*others[2]}
        eqs=[]
        for f in F:
            g=sp.Poly(sp.expand(f.subs(sub)),*others)
            eqs+=g.coeffs()
        sol=sp.solve(eqs,[n1,n2,n3],dict=True)
        sols.append((vars_[k],sol))
    print('s=',s)
    for v,sol in sols: print('  chart',v,'=lin(others):',sol)
