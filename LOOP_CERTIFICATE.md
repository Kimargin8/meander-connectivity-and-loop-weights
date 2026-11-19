# Reproducing the large-weight certificate

The theorem concerns ordered pairs of noncrossing matchings on \(2n\) labelled points and every fixed real loop weight \(q\ge18\). The proofs and normalization are in **large_loop_weight.pdf**. Independent repetition of the arithmetic certificate requires Python 3.10 or later and its standard library.

From the extracted directory run

~~~sh
python -B certificate18.py > replay.json
python -c "import json; assert json.load(open('replay.json')) == json.load(open('certificate18.json'))"
~~~

Run with assertions enabled. The program regenerates all 25,000 potential coordinates from the stated recurrence, checks the 128,600 extension inequalities, and checks the two state counts using a second combinatorial recurrence. A width-restricted forward count of returns gives a separate lower bound at the base weight. For small orders the full loop polynomials and the restrictions to irreducible and common-arch-free systems are also compared with direct matching pairs.

| File | Purpose |
|---|---|
| large_loop_weight.tex and windows.tex | Article and finite-window proof |
| certificate18.py | Exact upper and lower certificates and saddle inequalities |
| certificate18.json | Expected deterministic output |
| frontier.py | Canonical path configurations and geometric transitions |
| window_counts.py | Window counts by interval splitting and reflection |
| loop_coefficients.py | Finite loop polynomials and direct small diagrams |

The radius certificate is one-sided. The optimum over finite windows at fixed path weights controls uniform continuation growth from formal starting states. Equality with the empty-to-empty growth constant and the threshold of the large-weight phase remain open.
