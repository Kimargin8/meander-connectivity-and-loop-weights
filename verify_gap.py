"""Exact checks for the first distinct Smith gap. Python 3.10+, standard library."""
import argparse
from functools import cache
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
from math import comb, factorial, gcd
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent


def v2(x):
    assert x > 0
    return (x & -x).bit_length()-1


def pi_order(weight, degree):
    if not weight:
        return float("inf")
    return degree + 2*(v2(abs(weight.numerator))-v2(weight.denominator))


def top(k, h):
    n = 2*k+1
    assert n % h == 0
    m = n//h
    a = v2(k*(k+1))
    e = v2(factorial(k+1))
    triples = (comb(n, 3) + (2*n//3 if h%3 == 0 else 0))//h
    largest_count = 1 if k%2 else m+1
    next_count = m-1 if k <= 2 else (m*(k+1)-2 if k%2 else triples+m*(k-2)-2)
    return dict(k=k, N=n, h=h, top_exponent=e, top_count=largest_count,
                gap=a, next_exponent=e-a, next_count=next_count)


def check_histogram(histogram, k, h):
    case = top(k, h)
    assert sum(histogram) == comb(2*k+1, k)//h
    assert len(histogram)-1 == case["top_exponent"] and histogram[-1] == case["top_count"]
    b = case["next_exponent"]
    assert not any(histogram[b+1:-1]) and histogram[b] == case["next_count"]


@cache
def interior(points):
    if not points:
        return ((),)
    return tuple(((points[0],points[j]),)+a+b
                 for j in range(1,len(points),2)
                 for a in interior(points[1:j]) for b in interior(points[j+1:]))


def is_connected(left,right,n):
    parent = list(range(n))
    def root(a):
        while a != parent[a]:
            a = parent[a]
        return a
    for a,b in left+right:
        a,b = root(a),root(b)
        if a == b:
            return False
        parent[a] = b
    return True


def spin_terms(diagram):
    terms={}
    for choices in product(range(2),repeat=len(diagram)):
        mask=sum(1<<edge[side] for edge,side in zip(diagram,choices))
        # Coefficient for a left endpoint is -i, for a right endpoint 1.
        terms[mask]=(sum(choices)-len(diagram))%4
    return terms


def spin_check(k):
    n=2*k+1
    shapes=interior(tuple(range(1,n)))
    checks=0
    for q in shapes:
        for z in range(n):
            other=tuple(((a+z)%n,(b+z)%n) for a,b in q)
            right=spin_terms(other)
            for p in shapes:
                phase=[0]*4
                for mask,exponent in spin_terms(p).items():
                    if mask in right:
                        phase[(k+exponent+right[mask])%4]+=1
                assert phase[1]==phase[3]
                assert phase[0]-phase[2]==is_connected(p,other,n)
                checks+=1
    return checks


def unit_elimination(A, bits):
    """General row and column operations over Z/(2^bits)."""
    counts = [0]
    while A:
        mask = (1<<bits)-1
        pos = next(((i,j) for i,row in enumerate(A) for j,x in enumerate(row) if x%2),None)
        if pos is None:
            counts.append(0)
            A = [[x//2 for x in row] for row in A]
            bits -= 1
            assert bits > 0
            continue
        i,j = pos
        row = A.pop(i)
        inv = pow(row.pop(j),-1,mask+1)
        for a in A:
            c = a.pop(j)*inv & mask
            if c:
                for col,x in enumerate(row):
                    a[col] = (a[col]-c*x) & mask
        counts[-1] += 1
    return counts


def direct(k):
    n = 2*k+1
    shapes = interior(tuple(range(1,n)))
    phases = [[[int(is_connected(p,tuple(((a+z)%n,(b+z)%n) for a,b in q),n))
                for z in range(n)] for q in shapes] for p in shapes]
    rows=[]
    for h in range(1,n+1):
        if n%h:
            continue
        m=n//h
        A=[[sum(phases[i][j][(b-a+z)%n] for z in range(0,n,m))
            for j in range(len(shapes)) for b in range(m)]
           for i in range(len(shapes)) for a in range(m)]
        hist=unit_elimination(A,top(k,h)["top_exponent"]+2)
        check_histogram(hist,k,h)
        rows.append(dict(k=k,h=h,geometric_histogram=hist))
    return rows


def exact_rotation_orders(k):
    """Project labelled geometric matrices, filling other rotation orders by I."""
    n = 2*k+1
    e = top(k,1)["top_exponent"]
    gap = v2(k*(k+1))
    bits = e+2
    modulus = 1<<bits
    shapes = interior(tuple(range(1,n)))
    shifts = [[tuple(((a+z)%n,(b+z)%n) for a,b in q) for z in range(n)] for q in shapes]
    phases = [[[int(is_connected(p,rotated,n)) for rotated in family]
               for family in shifts] for p in shapes]
    divisors = [r for r in range(1,n+1) if n%r==0]
    projectors = {}
    triple_orbits = (comb(n,3)-(n//3 if n%3==0 else 0))//n
    result = []
    for r in divisors:
        # Start with averaging over <R^r>, then remove every proper-order part.
        coefficient = pow(n//r,-1,modulus)
        proj = [(coefficient if z%r==0 else 0) for z in range(n)]
        for d, earlier in projectors.items():
            if r%d==0:
                proj = [(x-y)%modulus for x,y in zip(proj,earlier)]
        projectors[r] = proj
        phi = sum(gcd(z,r)==1 for z in range(1,r+1))
        delta = int(r==1)
        short = int(n%3==0 and (n//3)%r==0)
        top_count = delta if k%2 else phi+delta
        next_count = ((k+1 if k%2 else triple_orbits+k-2+short)*phi-2*delta)
        # This is AP+I-P on the full labelled space. The complement adds only units.
        blocks = []
        for i,row in enumerate(phases):
            block_row = []
            for j,phase in enumerate(row):
                projected = [sum(proj[z]*phase[(b+z)%n] for z in range(n)) for b in range(n)]
                if i==j:
                    projected = [x+int(b==0)-proj[b] for b,x in enumerate(projected)]
                block_row.append([x%modulus for x in projected])
            blocks.append(block_row)
        matrix = [[blocks[i][j][(b-a)%n] for j in range(len(shapes)) for b in range(n)]
                  for i in range(len(shapes)) for a in range(n)]
        hist = unit_elimination(matrix,bits)
        hist.extend([0]*(e+1-len(hist)))
        assert hist[e]==top_count and hist[e-gap]==next_count
        assert not any(hist[e-gap+1:e])
        result.append(dict(k=k,rotation_order=r,top_count=top_count,next_count=next_count,
                           measured_positive_histogram=hist[1:]))
    return result


def avoids(k,u):
    return Fraction(comb(k+1,u),comb(2*k+1,u)) if u<=k+1 else Fraction(0)


def scalar_checks(limit):
    total=0
    for k in range(3,limit+1):
        a=v2(k*(k+1))
        for u in range(1 if k%2 else 2,2*k+2):
            B=avoids(k,u)
            depth=pi_order(B if k%2 else B-avoids(k,2),u)
            assert depth>=2*a
            assert (depth==2*a)==(u in ((2,4) if k%2 else (4,6)))
            total+=1
    return total


def residues(k,s,t):
    a=v2(k*(k+1))
    f=[]
    for j in range(min(s,t)+1):
        u=s+t-j
        B=avoids(k,u)
        error=B-avoids(k,s)*avoids(k,t) if k%2 else B-avoids(k,2)
        depth=pi_order(error,s+t)
        assert depth>=2*a
        f.append(int(depth==2*a))
    # Expand f(intersection size) in the binomial basis over F2.
    return tuple(j for j in range(len(f))
                 if sum(comb(j,l)*f[l] for l in range(j+1))%2)


def orbit_data(n,h,size,needed):
    m=n//h
    limit=(1<<n)-1
    reps=[]
    seen=set()
    incidence={j:{} for j in needed}
    for labels in combinations(range(n),size):
        mask=sum(1<<x for x in labels)
        if mask in seen:
            continue
        column=1<<len(reps)
        reps.append(labels)
        orbit={((mask<<z)&limit)|(mask>>(n-z)) for z in range(0,n,m)}
        seen.update(orbit)
        for member in orbit:
            points=[]
            while member:
                bit=member & -member
                points.append(bit.bit_length()-1)
                member^=bit
            for j,table in incidence.items():
                for sub in combinations(points,j):
                    table[sub]=table.get(sub,0)^column
    return reps,incidence


def binary_rank(rows):
    basis={}
    for row in rows:
        while row:
            lead=row.bit_length()-1
            if lead not in basis:
                basis[lead]=row
                break
            row^=basis[lead]
    return len(basis)


def incidence_check(k,h):
    n=2*k+1
    sizes=(1,2,3) if k%2 else (2,3,4)
    coeff={(s,t):residues(k,s,t) for s in sizes for t in sizes}
    data={t:orbit_data(n,h,t,{j for s in sizes for j in coeff[s,t]}) for t in sizes}
    offset={}
    width=0
    for t in sizes:
        offset[t]=width
        width+=len(data[t][0])
    rows=[]
    for s in sizes:
        for labels in data[s][0]:
            row=0
            for t in sizes:
                part=0
                for j in coeff[s,t]:
                    for sub in combinations(labels,j):
                        part^=data[t][1][j].get(sub,0)
                row |= part<<offset[t]
            rows.append(row)
    if k%2==0:
        # Restrict both arguments to differences from the first pair.
        p=len(data[2][0])
        ref=rows[0]
        rows=[(row^ref if i<p else row) for i,row in enumerate(rows) if i]
        rows=[(row^((1<<p)-1) if row&1 else row) & ~1 for row in rows]
    got=binary_rank(rows)
    case=top(k,h)
    assert got==case["next_count"]
    return dict(**case,generators=width,measured_binary_rank=got)


def previous_tables(directory):
    inputs=[directory/"cyclic_smith.json",directory/"subgroup_geometry.json"]
    count=0
    for path in inputs:
        for case in json.loads(path.read_bytes())["cases"]:
            if not case["k"]:
                continue
            if "h" in case:
                check_histogram(case["valuation_counts"],case["k"],case["h"])
            else:
                check_histogram(case["valuations"],case["k"],case["sites"])
            count+=1
    return {"cases":count,"input_sha256":{p.name:sha256(p.read_bytes()).hexdigest() for p in inputs}}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    record=dict(exact_scalar_cases=scalar_checks(256),previous_geometry=previous_tables(HERE),
                spin_graph_comparisons=sum(spin_check(k) for k in range(1,5)),
                independent_geometry=[c for k in range(1,6) for c in direct(k)],
                exact_rotation_orders=[c for k in (3,4) for c in exact_rotation_orders(k)],
                leading_forms=[incidence_check(k,h) for k,h in
                               [(4,3),(7,3),(8,1),(15,31),(16,11),(31,63),(32,65)]])
    text=json.dumps(record,indent=2)+"\n"
    if args.output:
        args.output.write_text(text,encoding="utf8")
    else:
        print(text,end="")


if __name__=="__main__":
    main()
