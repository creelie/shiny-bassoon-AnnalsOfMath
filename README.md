# Torsion D-Brane Charges from H-Flux on Calabi–Yau Threefolds

LaTeX source and verification code for the paper
*Torsion D-Brane Charges from H-Flux on Calabi–Yau Threefolds*
by Deep Bhattacharjee and Ushashi Bhattacharya.

## Layout

| Path | Contents |
| --- | --- |
| `tex/main.tex` | the paper (amsart, pdflatex) |
| `main.pdf` | the compiled paper |
| `code/twisted_k.py` | exact computation of the twisted K-groups of the paper (sympy, Smith normal forms) |
| `code/results.json` | the groups written by the script |

## Running the checks

```
cd code
python3 twisted_k.py
```

The script computes, from the cohomology data alone, the untwisted K-groups of
Table 1, the twisted K-groups of the closed 3-manifolds and of the simply
connected 6-manifolds of Theorems 1.1 and 1.2 (quintic, bicubic and Schoen
threefolds among the examples), the T-duality exchange of Proposition 5.2 for
circle bundles over surfaces, and the Grothendieck–Riemann–Roch and
Ramond–Ramond integrality identities of Section 6. It ends with
`31 checks passed, 0 failed`. Requires Python 3 and sympy.

## Building the paper

```
cd tex
pdflatex main.tex && pdflatex main.tex
```

## Contact

Deep Bhattacharjee, egspl@proton.me
