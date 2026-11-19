"""A rational window majorant and a width-restricted return lower bound."""
from fractions import Fraction as F
from functools import cache
from hashlib import sha256
from itertools import product
import json
from math import lcm, prod

from frontier import canonical, truncate, successors
from loop_coefficients import direct_rows, finite_data
from window_counts import count_windows


def windows(depth):
    """Generate matching sections directly, forbidding a pair outside the window."""
    result = set()
    for n in range(2*depth+1):
        for cut in range(n+1):
            lo,hi = max(0,cut-depth),min(2*n,cut+depth)
            @cache
            def intervals(first,end):
                if first==end:
                    return ((),)
                ans = []
                for second in range(first+1,end,2):
                    if not (lo<=first<hi or lo<=second<hi):
                        continue
                    for inside in intervals(first+1,second):
                        for outside in intervals(second+1,end):
                            ans.append(((first,second),)+inside+outside)
                return tuple(ans)
            for pairs in intervals(0,2*n):
                tags = [0]*(2*n)
                for label,(a,b) in enumerate(pairs):
                    tags[a] = tags[b] = label
                choices = [(1,3) if a<cut<=b else (2,) for a,b in pairs]
                for sizes in product(*choices):
                    result.add(canonical(tags[:cut][::-1],tags[cut:],sizes))
            intervals.cache_clear()
    return sorted(result)


def upper():
    weights = [1,F(1432,1000),F(2940,1000),F(4284,1000)]
    @cache
    def weight(tag):
        close,insert,remove = tag
        return 18**close*prod((weights[i] for i in insert),start=F(1))/prod(weights[i] for i in remove)
    states = windows(6)
    ids = {state:i for i,state in enumerate(states)}
    rows = [[] for _ in states]
    extensions = windows(7)
    for ext in extensions:
        i = ids[truncate(ext,6)]
        rows[i].append([(ids[truncate(dst,6)],weight(tag)) for _,dst,tag in successors(ext)])
    tags = {v for variants in rows for row in variants for _,v in row}
    scale = lcm(*(x.denominator for x in tags))
    for i,variants in enumerate(rows):
        rows[i] = [[(j,int(v*scale)) for j,v in row] for row in variants]
    def apply(vector):
        return [max(sum(vector[j]*v for j,v in row) for row in variants) for variants in rows]
    unit = 10**10
    h = [unit]*len(states)
    for _ in range(180):
        v = apply(h)
        peak = max(v)
        h = [max(1,unit*x//peak) for x in v]
    worst = max(F(x,scale*y) for x,y in zip(apply(h),h))
    assert count_windows(6)["orbits"]==len(states)
    assert count_windows(7)["orbits"]==len(extensions)
    assert worst<F(5447,1000), float(worst)
    return dict(windows=len(states),extensions=len(extensions),largest=str(worst),
                weights=list(map(str,weights[1:])),potential_sha256=sha256(
                    json.dumps(h,separators=(",",":")).encode()).hexdigest())


def returns(max_ends,order):
    states = [((),(),())]
    ids = {states[0]:0}
    rows = []
    for s in states:
        row = []
        for _,t,tag in successors(s):
            if len(t[0])+len(t[1])>max_ends:
                continue
            if t not in ids:
                ids[t]=len(states)
                states.append(t)
            row.append((ids[t],18**tag[0]))
        rows.append(row)
    v = [1]+[0]*(len(states)-1)
    counts = [1]
    for step in range(1,2*order+1):
        after = [0]*len(v)
        for x,row in zip(v,rows):
            if x:
                for j,f in row:
                    after[j] += f*x
        v = after
        if step%2==0:
            counts.append(v[0])
    return counts,len(states)


def main():
    h = upper()
    t,b = F(6065,10000),F(6977,10000)
    assert t*F(5447,1000)**2<18
    counts,nstates = returns(14,40)
    full,irre,free = finite_data()
    for a in range(1,7):
        assert direct_rows(a)==(full[a],irre[a],free[a])
    for a in range(8):
        # A history of length 2a can use at most 2a ends before returning.
        assert counts[a]==sum(x*18**(a-s) for s,x in free[a].items())
    v = sum(x*(t/18)**a for a,x in enumerate(counts))
    direct = t*sum((2*a+1)*x*(t/18)**a for a,x in enumerate(counts))-1
    terms = [(a,s,(2*a-1)*b**a*F(x,18**s)) for a,row in enumerate(irre) if a>=2 for s,x in row.items()]
    derivative = b-sum((2*s-a)*v for a,s,v in terms)
    sign = b-1+sum(v for _,__,v in terms)
    assert t*v*v>b and sign>0 and derivative>0 and direct>0
    print(json.dumps(dict(base=18,design=h,active_end_cap=14,lower_states=nstates,
                         orders=40,t=str(t),b=str(b),truncated_image=str(t*v*v),
                         derivative_identity_margin=str(direct),ray_sign=str(sign),
                         ray_derivative=str(derivative),closed_weighted_counts=counts),indent=2))


if __name__=="__main__":
    main()
