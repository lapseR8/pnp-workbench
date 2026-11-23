from modlin import *
# cross r=0,c=0 : basis
def cross(r=0,c=0,variant='a'):
    Bs=[]
    others_c=[j for j in range(Q) if j!=c]; others_r=[i for i in range(P) if i!=r]
    if variant=='a':
        # row sum off c =0: t(r,j)-t(r,j0)
        j0=others_c[0]
        for j in others_c[1:]:
            v=[0]*N; v[idx(r,j)]=1; v[idx(r,j0)]=-1; Bs.append(v)
        for i in others_r:
            v=[0]*N; v[idx(i,c)]=1; v[idx(r,c)]=-1; Bs.append(v)
    else:
        i0=others_r[0]
        for i in others_r[1:]:
            v=[0]*N; v[idx(i,c)]=1; v[idx(i0,c)]=-1; Bs.append(v)
        for j in others_c:
            v=[0]*N; v[idx(r,j)]=1; v[idx(r,c)]=-1; Bs.append(v)
    return Bs
for var in ("ab" if __name__=="__main__" else ""):
    Bs=cross(variant=var)
    print(var,len(Bs),valid_mod(Bs),'tangent dim',tangent_dim(Bs))
