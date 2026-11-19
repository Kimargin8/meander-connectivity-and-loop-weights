"""Exact finite windows of decorated paths at a cut through two matchings.

The upper and lower tuples list ends from the river outwards. Identical
integers are ends of the same open path. Reflection exchanges the lists.
"""
from functools import cache
from itertools import product


def canonical(upper, lower, sizes, mirror=True):
    """Rename paths in encounter order; choose the smaller reflected state."""
    def oriented(u, v):
        ids = dict.fromkeys((*u, *v))
        labels = {key: i for i, key in enumerate(ids)}
        return (tuple(labels[x] for x in u), tuple(labels[x] for x in v),
                tuple(sizes[x] for x in ids))
    return min(oriented(upper, lower), oriented(lower, upper)) if mirror else oriented(upper,lower)


@cache
def pairs(n):
    if n == 0:
        return ((),)
    return tuple(((0, 2*j+1), *((a+1,b+1) for a,b in inside),
                  *((a+2*j+2,b+2*j+2) for a,b in outside))
                 for j in range(n) for inside in pairs(j) for outside in pairs(n-j-1))


def representatives(depth):
    """All path pairings supported by at least one end in the window."""
    result = set()
    for paths in range(2*depth+1):
        width = 2*paths
        for cut in range(paths+1):  # the other cuts follow by reflection
            left, right = max(0,cut-depth), min(width,cut+depth)
            for matching in pairs(paths):
                if any(not (left<=i<right or left<=j<right) for i,j in matching):
                    continue
                labels = [0]*width
                for idx, (i,j) in enumerate(matching):
                    labels[i] = labels[j] = idx
                parities = [(1,3) if i<cut<=j else (2,) for i,j in matching]
                for lengths in product(*parities):
                    result.add(canonical(labels[:cut][::-1], labels[cut:], lengths))
    return sorted(result)


def truncate(state, depth, mirror=True):
    u, v, lengths = state
    retained = set(u[:depth]) | set(v[:depth])
    return canonical([x for x in u if x in retained], [x for x in v if x in retained], lengths,mirror)


def successors(state,mirror=True):
    """A multiplier is given as (cycle count, created classes, deleted classes)."""
    u, v, sizes = state
    n = len(sizes)
    yield "U",canonical((n,*u), (n,*v), (*sizes,1),mirror), (0,(1,),())
    for op,take, give in (("R",u,v),("L",v,u)):
        if take:
            p = take[0]
            old = sizes[p]
            updated = list(sizes)
            updated[p] = 2+((old+1)%2)
            uu,vv=take[1:],(p,*give)
            if op=="L":
                uu,vv=vv,uu
            yield op,canonical(uu, vv, updated,mirror), (0,(updated[p],),(old,))
    if not u or not v:
        return
    p, q = u[0], v[0]
    if p == q:
        if sizes[p] == 3:
            yield "D",canonical(u[1:],v[1:],sizes,mirror), (1,(),(3,))
        return
    new = 2+(sizes[p]+sizes[q]+1)%2
    updated = list(sizes)
    updated[p] = new
    uu = tuple(p if x==q else x for x in u[1:])
    vv = tuple(p if x==q else x for x in v[1:])
    yield "D",canonical(uu,vv,updated,mirror), (0,(new,),tuple(sorted((sizes[p],sizes[q]))))


def graph(depth):
    states = representatives(depth)
    ids = {s:i for i,s in enumerate(states)}
    extensions = representatives(depth+1)
    rows = [set() for _ in states]
    for ext in extensions:
        s = ids[truncate(ext,depth)]
        row = tuple(sorted((ids[truncate(dst,depth)],weight)
                           for op,dst,weight in successors(ext)))
        rows[s].add(row)  # multiplicities of the two slides are retained
    assert all(rows)
    return states,[sorted(variants) for variants in rows],len(extensions)
