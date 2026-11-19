"""Exact total, irreducible, and common-arch-free loop polynomials."""
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import product
from math import comb


def physical_options(state):
    cut,mate,lens=state
    size=len(mate)
    def inc(j):
        vals=list(lens);vals[j]+=1;vals[mate[j]]+=1
        return tuple(vals)
    shift=lambda k: k+2 if k>=cut else k
    ins=tuple(map(shift,mate))
    yield "U", (cut+1,ins[:cut]+(cut+1,cut)+ins[cut:],lens[:cut]+(1,1)+lens[cut:]),0
    if cut:
        yield "R",(cut-1,mate,inc(cut-1)),0
    if cut<size:
        yield "L",(cut+1,mate,inc(cut)),0
    if not 0<cut<size:
        return
    i,j=cut-1,cut;a,b=mate[i],mate[j];closed=a==j
    m=list(mate);l=list(lens)
    if not closed:
        m[a],m[b]=b,a
        l[a]=l[b]=lens[i]+lens[j]+1
    del m[i:j+1];del l[i:j+1]
    dest=(cut-1,tuple(k-2 if k>j else k for k in m),tuple(l))
    yield "D",dest,lens[i]+1 if closed else 0


def full_row(n):
    """Exact frontier count, using paired ports alone. Defect marks are decoded last."""
    width=(comb(2*n,n)**2).bit_length()
    mask=(1<<width)-1
    states={(0,()):1}
    for left in reversed(range(2*n)):
        future=defaultdict(int)
        for (c,m),count in states.items():
            # Lengths can be placeholders since each closure gets the same loop marker.
            for op,(d,mm,ll),cycle in physical_options((c,m,(1,)*len(m))):
                if max(d,len(mm)-d)<=left:
                    future[d,mm]+=count<<width if cycle else count
        states=future
    value=states[0,()]
    row={n-k:(value>>(k*width))&mask for k in range(1,n+1)}
    assert sum(row.values())==(comb(2*n,n)//(n+1))**2
    return row


def power_coefficient(a,m,n):
    terms=[{0:1}]+[{} for _ in range(n)]
    for _ in range(m):
        new=[Counter() for _ in range(n+1)]
        for k,poly in enumerate(terms):
            for d,x in poly.items():
                for j in range(min(n-k,len(a)-1)+1):
                    for s,y in a[j].items():new[k+j][s+d]+=x*y
        terms=new
    return terms[n]


def finite_data():
    N=12
    full=[{0:1}]+[full_row(n) for n in range(1,N+1)]
    irre=[{0:1},{0:1}]
    for n in range(2,N+1):
        vals=Counter(full[n])
        for i in range(1,n):
            inner=power_coefficient(full,2*i,n-i)
            for s,x in irre[i].items():
                for d,y in inner.items():vals[s+d]-=x*y
        assert all(v>=0 for v in vals.values())
        irre.append({s:v for s,v in vals.items() if v})
        assert all(i>=n/2 for i in irre[n])
    nobubbles=[{0:1},{}]
    nj=[{0:1},{}]+irre[2:]
    for n in range(2,N+1):
        st=power_coefficient(nj,2*n+1,n)
        assert all(v%(2*n+1)==0 for v in st.values())
        nobubbles.append({d:v//(2*n+1) for d,v in st.items()})
    return full,irre,nobubbles


@lru_cache(None)
def arches(n):
    if not n:return ((),)
    return tuple((r+1,*(x+1 for x in inn),0,*(x+r+2 for x in out))
                 for k in range(n) for r in [2*k]
                 for inn in arches(k) for out in arches(n-k-1))


def direct_rows(n):
    full,bubble,irre=Counter(),Counter(),Counter()
    for top,bottom in product(arches(n),repeat=2):
        remaining=set(range(2*n));loops=0
        while remaining:
            v=next(iter(remaining)); loops+=1
            while v in remaining:
                remaining.remove(v);remaining.remove(top[v]);v=bottom[top[v]]
        d=n-loops;full[d]+=1
        if all(top[j]!=bottom[j] for j in range(2*n)):bubble[d]+=1
        if not any(all(i<=top[v]<j and i<=bottom[v]<j for v in range(i,j))
                   for i in range(2*n) for j in range(i+2,2*n+1,2)
                   if j-i<2*n):
            irre[d]+=1
    return full,irre,bubble


