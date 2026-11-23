# Generate Singular scripts: equivariant (cyclic prime-order g) search for valid L of given dim, random generic chart per dim-vector.
import itertools, sys, os, random
from sympy.polys.rings import ring
from equiv import *
def epolys(g,p,z,orbs):
    names=[f'y{oi}_{k}' for oi,o in enumerate(orbs) for k in (range(p) if len(o)==p else [0])]
    R,*gens=ring(','.join(names),F); G=dict(zip(names,gens)); yn={'ring':R}
    for oi,o in enumerate(orbs):
        for k in (range(p) if len(o)==p else [0]): yn[(oi,k)]=G[f'y{oi}_{k}']
    Es=E_in_y(g,p,z,orbs,yn)
    def ps(e):
        out=[]
        for mon,c in e.terms():
            c=int(c)%q
            if c: out.append(str(c)+''.join(f'*{names[i]}^{ex}' for i,ex in enumerate(mon) if ex))
        return '+'.join(out) if out else '0'
    return names,[ps(e) for e in Es]
def script(typ,total,dims,ci,g,p,z,orbs,Vk,names,Estr,fname):
    random.seed(hash((typ,total,dims,ci))%1000003)
    pn=[];xn=[]
    for k in range(p):
        n=len(Vk[k])
        for a in range(dims[k]):
            xn.append(f'x{k}_{a}')
            for l in range(dims[k],n): pn.append(f'p{k}_{a}_{l}')
    if not pn: pn=['dummy']
    lines=[f'ring Ry={q},({",".join(names)}),dp;']
    for m,e in enumerate(Estr): lines.append(f'poly E{m+1}={e};')
    lines.append(f'ring R={q},({",".join(xn)},{",".join(pn)}),dp;')
    img={}
    for k in range(p):
        n=len(Vk[k]); Gm=[[random.randrange(1,q) for j in range(n)] for l in range(n)]
        for j,o in enumerate(Vk[k]):
            oi=orbs.index(o); terms=[]
            for a in range(dims[k]):
                for l in range(n):
                    if l<dims[k]:
                        if l==a: terms.append(f'{Gm[l][j]}*x{k}_{a}')
                    else: terms.append(f'{Gm[l][j]}*x{k}_{a}*p{k}_{a}_{l}')
            img[f'y{oi}_{k}']='+'.join(terms) if terms else '0'
    lines.append('map phi=Ry,'+','.join(img[nm] for nm in names)+';')
    xprod='*'.join(xn)
    lines.append('ideal I; int i; matrix M; poly f;')
    for m in range(1,5):
        lines.append(f'f=phi(E{m}); if(f!=0){{M=coef(f,{xprod}); for(i=1;i<=ncols(M);i++){{I=I+M[2,i];}}}}')
    tag=os.path.basename(fname)
    lines.append(f'option(redSB); ideal G=slimgb(I); if(G[1]==1){{print("{tag} UNIT");}} else {{print("{tag} NONUNIT dim "+string(dim(std(G))));}}')
    lines.append('quit;')
    open(fname,'w').write('\n'.join(lines)+'\n')
if __name__=="__main__":
    typ=sys.argv[1]; total=int(sys.argv[2]); outdir=sys.argv[3]; nch=int(sys.argv[4]); os.makedirs(outdir,exist_ok=True)
    gs={'c5':(make_g([0,1,2,3],[1,2,3,4,0]),5,[1,2,3,4]),'r3c3':(make_g([1,2,0,3],[1,2,0,3,4]),3,[1,2]),
        'c3':(make_g([0,1,2,3],[1,2,0,3,4]),3,[1,2]),'r3':(make_g([1,2,0,3],[0,1,2,3,4]),3,[1,2])}
    g,p,units=gs[typ]; z,orbs,Vk=setup(g,p); names,Estr=epolys(g,p,z,orbs)
    seen=set(); n=0
    for dims in cases(p,Vk,total):
        canon=min(tuple(dims[(u*k)%p] for k in range(p)) for u in units)
        if canon in seen: continue
        seen.add(canon)
        for ci in range(nch):
            script(typ,total,dims,ci,g,p,z,orbs,Vk,names,Estr,f"{outdir}/{typ}_{total}_{''.join(map(str,dims))}_{ci}.sing"); n+=1
    print(typ,total,n,'files')
