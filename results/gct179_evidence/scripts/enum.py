import itertools
from collections import Counter
rows=range(4)
pairs=list(itertools.combinations(rows,2))
res=set()
for r in itertools.product(range(0,5),repeat=4):
  for u in itertools.product(range(0,4),repeat=4):
    if sum(u)>7 or any(u[i]>r[i] for i in rows): continue
    if any(7-r[i] < sum(u[j] for j in rows if j!=i) for i in rows): continue
    ranges=[]
    for (i,j) in pairs:
        k,l=[x for x in rows if x not in (i,j)]
        lo=max(r[i],r[j]); hi=min(7,r[i]+r[j],7-u[k]-u[l])
        ranges.append(range(lo,hi+1))
    for pv in itertools.product(*ranges):
        p={}
        for (i,j),v in zip(pairs,pv): p[(i,j)]=p[(j,i)]=v
        ok=True
        for i,j,k,l in itertools.permutations(rows):
            if 7+r[i]>p[(i,j)]+p[(i,k)]+u[l]: ok=False;break
        if ok:
            # canonical under S4 row perms
            best=None
            for s in itertools.permutations(rows):
                key=(tuple(r[s[i]] for i in rows),tuple(u[s[i]] for i in rows),tuple(p[(s[i],s[j])] for i,j in pairs))
                if best is None or key<best: best=key
            res.add(best)
print('patterns',len(res))
U=Counter(tuple(sorted(x[1],reverse=True)) for x in res)
print('u-multisets',sorted(U.items()))
R=Counter(tuple(sorted(x[0],reverse=True)) for x in res)
print('r-multisets',sorted(R.items()))
