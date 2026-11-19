"""Small exact integer Smith forms from polygon chords and alternating paths."""
from collections import Counter
from itertools import combinations
from math import comb, factorial
import json

import sympy
from sympy.matrices.normalforms import smith_normal_form


def diagrams(k):
    n = 2*k+1
    for mark in range(n):
        for opens in combinations(range(2*k), k):
            open_set = set(opens)
            partners = [-1]*n
            stack = []
            for i in range(2*k):
                site = (mark+i+1)%n
                if i in open_set:
                    stack.append(site)
                elif stack:
                    other = stack.pop()
                    partners[other], partners[site] = site, other
                else:
                    break
            else:
                assert not stack
                yield tuple(partners)


def rotate(row, shift):
    n = len(row)
    return tuple(-1 if row[(i-shift)%n]<0 else (row[(i-shift)%n]+shift)%n
                 for i in range(n))


def connected(row, col):
    site = row.index(-1)
    edges = (col, row)
    length = 1
    while True:
        site = edges[(length-1)%2][site]
        if site<0:
            return length == len(row)
        length += 1
        assert length<=len(row)


def verify(k, h):
    n = 2*k+1
    m = n//h
    states = set(diagrams(k))
    assert len(states) == comb(n,k)
    orbits = []
    while states:
        first = min(states)
        orbit = [rotate(first,m*j) for j in range(h)]
        assert len(set(orbit)) == h
        states.difference_update(orbit)
        orbits.append(orbit)
    A = sympy.Matrix([[sum(connected(row[0],col) for col in orbit) for orbit in orbits]
                      for row in orbits])
    assert A == A.T
    factors = Counter(abs(int(x)) for x in smith_normal_form(A,domain=sympy.ZZ).diagonal())
    assert 0 not in factors and all(x & (x-1) == 0 for x in factors)
    v2 = lambda x: (x & -x).bit_length()-1
    e = v2(factorial(k+1))
    gap = v2(k*(k+1))
    top = 2**e
    next_factor = top//(2**gap)
    triple_orbits = (comb(n,3)+(2*n//3 if h%3==0 else 0))//h
    expected_next = m-1 if k<=2 else (m*(k+1)-2 if k%2 else triple_orbits+m*(k-2)-2)
    assert max(factors)==top and factors[top] == (1 if k%2 else m+1)
    assert factors[next_factor]==expected_next
    assert not any(next_factor<x<top for x in factors)
    return dict(k=k,subgroup_order=h,dimension=A.rows,integer_factors=dict(sorted(factors.items())))


if __name__ == "__main__":
    cases = [verify(k,h) for k in range(1,6) for h in range(1,2*k+2)
             if (2*k+1)%h==0 and comb(2*k+1,k)//h<=50]
    print(json.dumps(dict(sympy=sympy.__version__,cases=cases),indent=2))
