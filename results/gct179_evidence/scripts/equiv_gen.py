# generate Singular files for equivariant cases (same math as equiv.py)
import itertools, sys, os
import equiv
from equiv import *
def gen_case(g,p,dims,charts,z,orbs,Vk,fname,GR):
    pnames=[];xinames=[]
    for k in range(p):
        d=dims[k]; n=len(Vk[k])
        for a in range(d):
            xinames.append(f'x{k}_{a}')
            for j in range(n):
                if j not in charts[k]: pnames.append(f'p{k}_{a}_{j}')
    R,*gens=ring(','.join(pnames+xinames),F)
    G=dict(zip(pnames+xinames,gens)); yn={'ring':R}
    for k in range(p):
        for j,o in enumerate(Vk[k]):
            oi=orbs.index(o); val=R.zero
            for a in range(dims[k]):
                piv=charts[k]
                for l in range(len(Vk[k])):
                    coef = R.one if (l in piv and piv.index(l)==a) else (R.zero if l in piv else G[f'p{k}_{a}_{l}'])
                    if coef!=R.zero and GR[k][l][j]: val+=coef*G[f'x{k}_{a}']*GR[k][l][j]
            yn[(oi,k)]=val
    for oi,o in enumerate(orbs):
        for k in range(p):
            if (oi,k) not in yn and (len(o)==p or k==0): yn[(oi,k)]=R.zero
    Es=E_in_y(g,p,z,orbs,yn); npar=len(pnames)
    eqs=set()
    for e in Es:
        coll={}
        for mon,c in e.terms():
            coll.setdefault(mon[npar:],{}); d=coll[mon[npar:]]; d[mon[:npar]]=d.get(mon[:npar],0)+int(c)
        for xm,d in coll.items():
            poly=tuple(sorted((pm,c%q) for pm,c in d.items() if c%q))
            if poly: eqs.add(poly)
    if any(all(sum(pm)==0 for pm,c in poly) for poly in eqs): return 'unit(const)'
    if npar==0: return 'INCONSISTENT' if eqs else 'VALID'
    def ps(poly):
        return '+'.join(f"{c}"+''.join(f"*{pnames[i]}^{e}" for i,e in enumerate(pm) if e) for pm,c in poly)
    with open(fname,'w') as f:
        f.write(f"ring r={q},({','.join(pnames)}),dp;\noption(redSB);\nideal I=\n"+",\n".join(ps(pl) for pl in eqs)+";\n")
        f.write(f'ideal G=slimgb(I);\nif(G[1]==1){{print("{os.path.basename(fname)} UNIT");}}else{{print("{os.path.basename(fname)} NONUNIT dim "+string(dim(std(G))));}}\nquit;\n')
    return 'FILE'
if __name__=="__main__":
    typ=sys.argv[1]; total=int(sys.argv[2]); outdir=sys.argv[3]; os.makedirs(outdir,exist_ok=True)
    gs={'c5':(make_g([0,1,2,3],[1,2,3,4,0]),5,[1,2,3,4]),'r3c3':(make_g([1,2,0,3],[1,2,0,3,4]),3,[1,2]),
        'c3':(make_g([0,1,2,3],[1,2,0,3,4]),3,[1,2]),'r3':(make_g([1,2,0,3],[0,1,2,3,4]),3,[1,2])}
    g,p,units=gs[typ]; z,orbs,Vk=setup(g,p)
    seen=set(); stats={}
    for dims in cases(p,Vk,total):
        canon=min(tuple(dims[(u*k)%p] for k in range(p)) for u in units)
        if canon in seen: continue
        seen.add(canon)
        import random
        nch=int(sys.argv[4]) if len(sys.argv)>4 else 1
        for ci in range(nch):
            random.seed(hash((typ,total,dims,ci))%100000)
            GR={k:[[random.randrange(q) for j in range(len(Vk[k]))] for l in range(len(Vk[k]))] for k in range(p)}
            charts=[tuple(range(dims[k])) for k in range(p)]
            fn=f"{outdir}/{typ}_{total}_{''.join(map(str,dims))}_{ci}.sing"
            r=gen_case(g,p,dims,charts,z,orbs,Vk,fn,GR)
            stats[r]=stats.get(r,0)+1
            if r=='VALID': print('VALID',dims,charts)
    print(typ,total,'dimV',[len(Vk[k]) for k in range(p)],stats)
