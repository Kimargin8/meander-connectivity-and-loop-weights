"""Check the integral-image step directly, using only labelled diagrams."""
from functools import cache
from itertools import combinations
from math import comb
import json


@cache
def matchings(sites):
    if not sites:
        return ((),)
    return tuple(((sites[0],sites[j]),)+a+b for j in range(1,len(sites),2)
                 for a in matchings(sites[1:j]) for b in matchings(sites[j+1:]))


def image(k):
    n = 2*k+1
    good = {}
    for size in range(k+1):
        for s in combinations(range(n),size):
            if all(x+1>=2*(i+1) for i,x in enumerate(s)):
                good[sum(1<<i for i in s)] = len(good)
    pivots = {}
    for mark in range(n):
        cut = tuple((mark+j)%n for j in range(1,n))
        for diagram in matchings(cut):
            subsets = [0]
            for a,b in diagram:
                # Each subset chooses zero or one end of a path.
                subsets += [s|(1<<end) for s in subsets for end in (a,b)]
            bits = sum(1<<good[s] for s in subsets if s in good)
            while bits:
                col = bits.bit_length()-1
                if col not in pivots:
                    pivots[col] = bits
                    break
                bits ^= pivots[col]
            assert bits
    assert len(pivots)==len(good)==comb(n,k)
    return dict(k=k,sites=n,binary_rank=len(pivots))


if __name__=="__main__":
    print(json.dumps([image(k) for k in range(7)],indent=2))
