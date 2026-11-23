# Exact (mod q) search for 7-dim valid L invariant under a prime-order cell permutation g.
import itertools, sys, time, signal
import sympy as sp
from sympy.polys.rings import ring
from sympy.polys.domains import GF
from core import MON, P, Q, N, idx
q=30091
F=GF(q)
def prim_root_of_unity(p):
    for a in range(2,q):
        z=pow(a,(q-1)//p,q)
        if z!=1: return z
def orbits(g):
    seen=set(); orbs=[]
    for v in range(N):
        if v in seen: continue
        o=[v]; w=g[v]
        while w!=v: o.append(w); w=g[w]
        seen.update(o); orbs.append(o)
    return orbs
def make_g(rowperm,colperm):
    return [idx(rowperm[v//Q],colperm[v%Q]) for v in range(N)]
def setup(g,p):
    z=prim_root_of_unity(p)
    orbs=orbits(g)
    # V_k basis: orbits of size p give one vector per k; fixed points only k=0
    Vk={k:[o for o in orbs if len(o)==p or k==0] for k in range(p)}
    return z,orbs,Vk
def E_in_y(g,p,z,orbs,ynames):
    # ynames[(orbit index,k)] -> generator; returns list of E_m polys in ring
    R=ynames['ring']
    tcell=[None]*N
    for oi,o in enumerate(orbs):
        for l,v in enumerate(o):
            if len(o)==1: tcell[v]=ynames[(oi,0)]
            else: tcell[v]=sum((ynames[(oi,k)]*pow(z,(-k*l)%p,q) for k in range(p)),R.zero)
    Es=[]
    for m in range(1,5):
        e=R.zero
        for mon in MON[m]:
            pr=R.one
            for v in mon: pr=pr*tcell[v]
            e+=pr
        Es.append(e)
    return Es
def run_case(g,p,dims,charts,z,orbs,Vk,timeout=60):
    # params for each k: chart pivots
    pnames=[];xinames=[]
    for k in range(p):
        d=dims[k]; n=len(Vk[k])
        for a in range(d):
            xinames.append(f'x{k}_{a}')
            for j in range(n):
                if j not in charts[k]: pnames.append(f'p{k}_{a}_{j}')
    R,*gens=ring(','.join(pnames+xinames),F)
    G=dict(zip(pnames+xinames,gens))
    yn={'ring':R}
    for k in range(p):
        for j,o in enumerate(Vk[k]):
            oi=orbs.index(o)
            val=R.zero
            for a in range(dims[k]):
                piv=charts[k]
                coef = R.one if (j in piv and piv.index(j)==a) else (R.zero if j in piv else G[f'p{k}_{a}_{j}'])
                val+=coef*G[f'x{k}_{a}']
            yn[(oi,k)]=val
    for oi,o in enumerate(orbs):
        for k in range(p):
            if (oi,k) not in yn and (len(o)==p or k==0): yn[(oi,k)]=R.zero
    Es=E_in_y(g,p,z,orbs,yn)
    nx=len(xinames); npar=len(pnames)
    eqs=set()
    for e in Es:
        coll={}
        for mon,c in e.terms():
            xm=mon[npar:]; pm=mon[:npar]
            coll.setdefault(xm,{}); coll[xm][pm]=(coll[xm].get(pm,0)+c)
        for xm,d in coll.items():
            poly=tuple(sorted((pm,int(c)%q) for pm,c in d.items() if int(c)%q))
            if poly: eqs.add(poly)
    if not npar:
        return ('const', 'INCONSISTENT' if eqs else 'VALID')
    pv=sp.symbols(pnames)
    exprs=[]
    for poly in eqs:
        ex=0
        for pm,c in poly: ex+=c*sp.prod([pv[i]**e for i,e in enumerate(pm)])
        exprs.append(ex)
    # constant equation?
    if any(all(sum(pm)==0 for pm,c in poly) for poly in eqs): return (npar,'unit(const)')
    def handler(s,f): raise TimeoutError
    signal.signal(signal.SIGALRM,handler); signal.alarm(timeout)
    try:
        GB=sp.groebner(exprs,*pv,modulus=q,order='grevlex')
        signal.alarm(0)
        return (npar, 'UNIT' if list(GB.exprs)==[1] else f'NONUNIT {len(GB.exprs)}')
    except TimeoutError:
        return (npar,'TIMEOUT')
def cases(p,Vk,total=7,maxd0=None):
    dimsV=[len(Vk[k]) for k in range(p)]
    for dims in itertools.product(*[range(min(dimsV[k],total)+1) for k in range(p)]):
        if sum(dims)!=total: continue
        ok=all(dims[k]+dims[(-k)%p]<=dimsV[k] for k in range(1,p))
        if maxd0 is not None and dims[0]>maxd0: ok=False
        if ok: yield dims
if __name__=="__main__":
    typ=sys.argv[1]
    if typ=='c5': g=make_g([0,1,2,3],[1,2,3,4,0]); p=5; units=[1,2,3,4]
    if typ=='r3c3': g=make_g([1,2,0,3],[1,2,0,3,4]); p=3; units=[1,2]
    if typ=='c3': g=make_g([0,1,2,3],[1,2,0,3,4]); p=3; units=[1,2]
    if typ=='r3': g=make_g([1,2,0,3],[0,1,2,3,4]); p=3; units=[1,2]
    total=int(sys.argv[2]) if len(sys.argv)>2 else 7
    z,orbs,Vk=setup(g,p)
    print(typ,'dim V_k',[len(Vk[k]) for k in range(p)],flush=True)
    seen=set()
    for dims in cases(p,Vk,total):
        # symmetry k->u*k (column/row perms normalizing <g>); take canonical rep
        canon=min(tuple(dims[(u*k)%p] for k in range(p)) for u in units)
        if canon in seen: continue
        seen.add(canon)
        res=[]
        chartlists=[list(itertools.combinations(range(len(Vk[k])),dims[k])) for k in range(p)]
        t0=time.time()
        for charts in itertools.product(*chartlists):
            r=run_case(g,p,dims,charts,z,orbs,Vk,timeout=int(sys.argv[3]) if len(sys.argv)>3 else 60)
            res.append(r[1])
        summary={s:res.count(s) for s in set(res)}
        print(dims,'charts',len(res),summary,'%.1fs'%(time.time()-t0),flush=True)
