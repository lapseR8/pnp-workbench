import numpy as np, sys
from tanext import T
from search import Prob
from invariants import inv
def Fb(b):
    a=b/(1-b); c=3*b-1
    V=[]
    def M(): return np.zeros((4,5),complex)
    x=M(); x[0,0]=1; x[1,0]=-1; V.append(x)
    x=M(); x[3,0]=1; x[1,0]=-1; V.append(x)
    x=M(); x[:,2]=[-1/3,-1/3,1,-1/3]; V.append(x)
    x=M(); x[0]=[0,1,0,-a,a-1]; V.append(x)
    x=M(); x[2]=[1,-b,0,b-1,0]; V.append(x)
    x=M(); x[2]=[0,-c,0,c-1,1]; V.append(x)
    return np.array([v.flatten() for v in V])
P=Prob(4,5,6)
for b in [0.2, 2/7, -3/11, 0.37+0.21j, 1.3-0.4j]:
    A=Fb(b)
    print('b=',b,'resid',np.linalg.norm(P.F(A)[0]),'dimT',T(A)[0], inv(A,4,5))
