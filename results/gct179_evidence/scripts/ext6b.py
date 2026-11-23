from common import *
exec(open('ext6.py').read().split("v=sp.symbols")[0])
vs=[-16,4,0,0,0, 0,1,1,1,1, 0,1,1,1,1, 0,1,1,1,1]
s=sp.Symbol('s')
tv=[t[c]+s*vs[c] for c in range(20)]
for k in range(1,5):
    print(k, sp.factor(Ek_sym(tv,k)))
