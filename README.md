# Meander connectivity and loop weights

Two preprints by Iaroslav Argunov, 19 November 2026. Each paper includes a proof and draws on a different model of meanders.

### Unbounded gaps between the two largest distinct Smith factors of a one-defect connectivity matrix

On \(2k+1\) labelled sites, mark one site and match the rest without crossings. For any subgroup of rotations, count the relative rotations under which two states glue to a connected path. For \(k\ge3\) the quotient of the two largest distinct Smith factors of this integer matrix is \(2^{v_2(k(k+1))}\). The paper determines both multiplicities. [Paper](gap.pdf) · [LaTeX](gap.tex) · [Rotation orders](rotation_orders.md)

### Loop counts in closed meandric systems at large weight

Give each loop in a pair of noncrossing perfect matchings weight \(q\). For every fixed real \(q\ge18\), the partition function has an \(n^{-3/2}\) prefactor. A finite window on the active ends gives a rational certificate for this interval. The paper also describes the Poisson approximation to the deficit in the number of loops as \(q\) increases. The location of the critical loop weight remains open. [Paper](large_loop_weight.pdf) · [LaTeX](large_loop_weight.tex) · [Window proof](windows.tex)

### Reproducing the computations

Extract the source archive or download the repository. Run these commands from its top directory using Python 3.10 or later, with assertions enabled. The reference environment used Python 3.13.5. All commands except the optional integer Smith calculation use only the standard library.

| Command | Output to compare |
|---|---|
| `python -B verify_gap.py --output replay.json` | `gap_checks.json` |
| `python -B image_check.py` | `image_check.json` |
| `python -B certificate18.py` | `certificate18.json` |
| `python -B integer_smith.py` | `integer_smith.json`, with SymPy 1.14.0 installed |

Compare parsed JSON objects. The optional integer Smith output records the version of SymPy used. The files `cyclic_smith.json` and `subgroup_geometry.json` are input profiles calculated from diagram gluings. The independent integer Smith check constructs the state matrices afresh using alternating paths and works over the integers. At larger orders the program checks the smaller incidence forms specified in the proof.

The loop certificate checks a positive vector on 25,000 windows against all 128,600 extensions, then counts width-restricted returns through order 40. The rational inequalities printed by the checker give the bounds used in the paper. The LaTeX sources build with Tectonic or a TeX distribution containing the AMS packages and Latin Modern.

### Citation and attribution

Cite the title of the paper you use, Iaroslav Argunov, the release date and the versioned release URL. Bibliographic entries are available in [CITATION.cff](CITATION.cff). The manuscripts are author preprints. AI tools assisted with research and the preparation and testing of verification programs. See [NOTICE.md](NOTICE.md) for sources and licenses.
