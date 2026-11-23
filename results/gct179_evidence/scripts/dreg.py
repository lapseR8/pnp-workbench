import sympy as sp
z=sp.symbols('z')
def dreg(n,degs,maxd=20):
    s=sp.series(sp.prod([(1-z**d) for d in degs])/(1-z)**n,z,0,maxd).removeO()
    c=[s.coeff(z,i) for i in range(maxd)]
    for i,x in enumerate(c):
        if x<=0: return i,[c[j] for j in range(i+1)]
for name,n,degs in [('3x5 L6 plain',48,[2]*21+[3]*56),('3x5 L6 Qadapt',27,[6]*56),('4x5 L7 plain',84,[2]*28+[3]*84+[4]*210),('4x5 Qadapt',56,[6]*84+[8]*210),('3x4 L5 plain',30,[2]*15+[3]*35),('3x3 L4 plain',16,[2]*10+[3]*20)]:
    d,c=dreg(n,degs); print(name,'dreg',d,'binom',sp.binomial(n+d,d))
