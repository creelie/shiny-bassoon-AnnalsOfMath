"""Exact computations for the paper
"Torsion D-Brane Charges from H-Flux on Calabi-Yau Threefolds".

Everything is computed from topological input (Betti numbers, torsion,
cup products) with exact integer and rational arithmetic; nothing is
copied from the text.  Requires sympy (for Smith normal forms); numpy and
matplotlib are not needed.

Run:  python3 twisted_k.py          (prints every result and writes results.json)
"""
from fractions import Fraction as F
from itertools import product
from math import comb, gcd
import json
import math

from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form

RESULTS = {}
PASSED = 0


def check(name, cond):
    global PASSED
    if not cond:
        raise AssertionError(name)
    PASSED += 1
    print(f"  ok  {name}")


# ----------------------------------------------------------------------
# 0. Finitely generated abelian groups and Smith normal form
# ----------------------------------------------------------------------
def group_from_relations(ngen, relations):
    """Invariants (free rank, [torsion orders]) of Z^ngen / <rows of relations>."""
    if ngen == 0:
        return (0, [])
    if not relations:
        return (ngen, [])
    M = Matrix(relations)
    S = smith_normal_form(M, domain=ZZ)
    diag = [abs(int(S[i, i])) for i in range(min(S.shape))]
    tors = [d for d in diag if d > 1]
    rank = ngen - sum(1 for d in diag if d != 0)
    return (rank, sorted(tors))


def kernel_rank_and_image(mat, nrows_src):
    """For an integer matrix A: Z^src -> Z^tgt (columns index the source),
    return the free rank of ker A and the rows spanning im A (as relations
    on the target)."""
    if nrows_src == 0:
        return 0, []
    A = Matrix(mat)
    r = A.rank()
    return A.cols - r, [list(A.col(j)) for j in range(A.cols)]


def fmt(g):
    rank, tors = g
    parts = []
    if rank == 1:
        parts.append("Z")
    elif rank > 1:
        parts.append(f"Z^{rank}")
    parts += [f"Z/{t}" for t in tors]
    return " + ".join(parts) if parts else "0"


# ----------------------------------------------------------------------
# 1. Untwisted K-theory from the Atiyah-Hirzebruch spectral sequence
#    for spaces whose integral cohomology is torsion-free and concentrated
#    in even degrees (every differential has source or target zero).
# ----------------------------------------------------------------------
def k_groups_even(betti):
    """betti: dict degree -> rank, all degrees even, torsion-free."""
    assert all(p % 2 == 0 for p in betti)
    return (sum(betti.values()), []), (0, [])


def k_groups_torus(n):
    """Kunneth: K^*(T^n) is the exterior algebra on n odd generators."""
    even = sum(comb(n, j) for j in range(0, n + 1, 2))
    odd = sum(comb(n, j) for j in range(1, n + 1, 2))
    return (even, []), (odd, [])


print("1. K-groups of standard spaces (computed)")
table = {}
table["S^2"] = k_groups_even({0: 1, 2: 1})
table["S^4"] = k_groups_even({0: 1, 4: 1})
table["CP^2"] = k_groups_even({0: 1, 2: 1, 4: 1})
table["CP^3"] = k_groups_even({0: 1, 2: 1, 4: 1, 6: 1})
table["K3"] = k_groups_even({0: 1, 2: 22, 4: 1})
table["T^2"] = k_groups_torus(2)
table["T^6"] = k_groups_torus(6)
for name, (k0, k1) in table.items():
    print(f"  {name:<6} K^0 = {fmt(k0):<8} K^1 = {fmt(k1)}")
check("S^2: K^0=Z^2, K^1=0", table["S^2"] == ((2, []), (0, [])))
check("S^4: K^0=Z^2, K^1=0", table["S^4"] == ((2, []), (0, [])))
check("K3: K^0=Z^24, K^1=0", table["K3"] == ((24, []), (0, [])))
check("T^6: K^0=Z^32, K^1=Z^32", table["T^6"] == ((32, []), (32, [])))
RESULTS["untwisted_table"] = {k: [fmt(a), fmt(b)] for k, (a, b) in table.items()}

# Bott periodicity: reduced K-theory of spheres
print("\n1b. Reduced K-theory of spheres from Bott periodicity")
for n in range(0, 8):
    k0t = 1 if n % 2 == 0 else 0
    k1t = 1 - k0t
    print(f"  S^{n}: K~^0 = {'Z' if k0t else 0}, K~^1 = {'Z' if k1t else 0}")

# ----------------------------------------------------------------------
# 2. Chern character and Todd class (symbolic, in Q[H]/(H^{n+1}))
# ----------------------------------------------------------------------
def poly_mul(a, b, n):
    return tuple(sum(a[i] * b[k - i] for i in range(k + 1) if i < len(a) and k - i < len(b))
                 for k in range(n + 1))


def poly_inv(a, n):
    """Inverse of a power series with constant term 1, modulo degree n+1."""
    inv = [F(1)] + [F(0)] * n
    for k in range(1, n + 1):
        inv[k] = -sum(a[i] * inv[k - i] for i in range(1, k + 1) if i < len(a))
    return tuple(inv)


def exp_series(c, n):
    """exp(c H) modulo H^{n+1}."""
    out, term = [], F(1)
    for k in range(n + 1):
        out.append(term)
        term = term * c / (k + 1)
    return tuple(out)


def todd_line(c, n):
    """Td of a line bundle with c_1 = c H:  cH/(1-e^{-cH})."""
    # Bernoulli expansion x/(1-e^{-x}) = 1 + x/2 + x^2/12 - x^4/720 + ...
    coeff = [F(1), F(1, 2), F(1, 12), F(0), F(-1, 720), F(0), F(1, 30240)]
    return tuple(coeff[k] * F(c) ** k if k < len(coeff) else F(0) for k in range(n + 1))


print("\n2. Chern character of a rank-4 bundle on CP^2 with c_1 = 2H, c_2 = H^2")
r, c1, c2 = 4, 2, 1
ch = (F(r), F(c1), F(c1 * c1 - 2 * c2, 2))
print(f"  ch = {ch[0]} + {ch[1]} H + {ch[2]} H^2")
check("ch(E) = 4 + 2H + H^2", ch == (4, 2, 1))

# ----------------------------------------------------------------------
# 3. Grothendieck-Riemann-Roch for the hyperplane CP^1 in CP^2 and for
#    the identity embedding; comparison of the two candidate formulas.
# ----------------------------------------------------------------------
print("\n3. Riemann-Roch for embeddings")
n = 2
# Koszul: i_! 1 = 1 - [O(-1)]  in K^0(CP^2),  ch = 1 - e^{-H}
koszul = tuple(F(int(k == 0)) - exp_series(-1, n)[k] for k in range(n + 1))
# i_*(1) = H,  i_*(h) = H^2 (h the generator of H^2(CP^1)).
def push(coeffs):
    return (F(0), coeffs[0], coeffs[1])
td_cp1 = todd_line(2, 1)            # Td(T CP^1) = 1 + h   (c_1 = 2h)
td_normal = todd_line(1, 1)         # N = O(1),  Td = 1 + h/2
correct = push(poly_inv(td_normal, 1))                   # i_*(ch(E)/Td(N))
wrong = push(poly_mul(td_cp1, poly_inv(td_normal, 1), 1))  # i_*(ch(E) Td(W)/Td(N))
print(f"  Koszul:            ch(i_! 1) = {koszul}")
print(f"  i_*(1/Td(N)):                  {correct}")
print(f"  i_*(Td(W)/Td(N)):              {wrong}")
check("GRR: ch(i_!1) = i_*(1/Td(N)) on CP^1 in CP^2", koszul == correct)
check("the formula with Td(W)/Td(N) fails on CP^1 in CP^2", koszul != wrong)
# Identity embedding of CP^1: N = 0, so the correct formula gives ch(E)=1,
# while Td(W)/Td(N) = Td(CP^1) = 1 + h.
check("identity embedding: Td(W)/Td(N) = 1 + h differs from 1", td_cp1 != (F(1), F(0)))
# Global form ch(i_! E) Td(X) = i_*(ch(E) Td(W)):
td_cp2 = poly_mul(todd_line(1, 2), poly_mul(todd_line(1, 2), todd_line(1, 2), 2), 2)  # (1+H/2+H^2/12)^3
lhs = poly_mul(koszul, td_cp2, 2)
rhs = push(td_cp1)
check("GRR global form ch(i_!1) Td(CP^2) = i_*(Td(CP^1))", lhs == rhs)
RESULTS["grr_hyperplane"] = {"koszul": [str(x) for x in koszul],
                             "i_*(1/Td N)": [str(x) for x in correct],
                             "i_*(Td W/Td N)": [str(x) for x in wrong]}

# Minasian-Moore / Freed-Hopkins: G = ch(x) sqrt(A-hat(X)) is not integral.
print("\n4. Non-integrality of the Ramond-Ramond field of O(1) on CP^2")
ch_line = exp_series(1, 2)
# A-hat(CP^2) = 1 - p_1/24,  p_1 = 3H^2  ->  sqrt = 1 - H^2/16
sqrt_ahat = (F(1), F(0), F(-1, 16))
G = poly_mul(ch_line, sqrt_ahat, 2)
print(f"  G(O(1)) = {G}")
check("degree-4 component of G(O(1)) is 7/16, not an integer", G[2] == F(7, 16))
RESULTS["rr_field_O1_CP2"] = [str(x) for x in G]

# ----------------------------------------------------------------------
# 5. Twisted K-theory from the twisted AHSS with d_3 = Sq^3_Z - H cup (-)
# ----------------------------------------------------------------------
print("\n5. Twisted K-theory of closed oriented 3-manifolds")


def twisted_k_3manifold(b1, tors_H2, k):
    """Closed oriented 3-manifold with H^1 = Z^b1, H^2 = Z^b1 + tors_H2,
    H^3 = Z and twist k times the generator of H^3.  Sq^3_Z vanishes on
    H^0, so d_3: H^0 -> H^3 is multiplication by k; all other d_3 vanish
    for degree reasons."""
    e0 = group_from_relations(1, [[k]]) if k == 0 else (0, [])   # ker(k)
    e2 = (b1, sorted(t for t in tors_H2 if t > 1))
    e1 = (b1, [])
    e3 = group_from_relations(1, [[k]])                          # Z/k
    # K^0 = E^0 (+) E^2 (E^0 free, so the extension splits);
    # K^1: 0 -> E^3 -> K^1 -> E^1 -> 0 with E^1 free, splits.
    k0 = (e0[0] + e2[0], e2[1])
    k1 = (e1[0] + e3[0], e3[1])
    return k0, k1


for name, b1, t, k in [("S^3", 0, [], 0), ("S^3", 0, [], 1), ("S^3", 0, [], 5),
                        ("T^3", 3, [], 0), ("T^3", 3, [], 4),
                        ("RP^3", 0, [2], 0), ("RP^3", 0, [2], 3),
                        ("L(6,1)", 0, [6], 0), ("L(6,1)", 0, [6], 4)]:
    k0, k1 = twisted_k_3manifold(b1, t, k)
    print(f"  {name:<7} k={k}:  K^0 = {fmt(k0):<12} K^1 = {fmt(k1)}")
    RESULTS.setdefault("three_manifolds", []).append([name, k, fmt(k0), fmt(k1)])
check("S^3, twist k: K^0 = 0, K^1 = Z/k", twisted_k_3manifold(0, [], 5) == ((0, []), (0, [5])))
check("S^3 untwisted: K^0 = Z, K^1 = Z", twisted_k_3manifold(0, [], 0) == ((1, []), (1, [])))
check("T^3, twist 4: K^0 = Z^3, K^1 = Z^3 + Z/4", twisted_k_3manifold(3, [], 4) == ((3, []), (3, [4])))

# Circle bundles M_e over a closed surface of genus g, Euler number e != 0:
# H^1 = Z^{2g},  H^2 = Z^{2g} + Z/e,  H^3 = Z  (Gysin sequence).
print("\n5b. T-duality check for circle bundles over surfaces (Bouwknegt-Evslin-Mathai)")


def circle_bundle_groups(g, e, k):
    if e == 0:
        return twisted_k_3manifold(2 * g + 1, [], k)
    return twisted_k_3manifold(2 * g, [abs(e)], k)


ok = True
for g in range(0, 3):
    for e in range(0, 6):
        for k in range(0, 6):
            k0, k1 = circle_bundle_groups(g, e, k)
            d0, d1 = circle_bundle_groups(g, k, e)    # T-dual: (e,k) -> (k,e)
            ok = ok and (k0 == d1) and (k1 == d0)
            if g == 1 and e in (0, 2) and k in (0, 3):
                print(f"  g={g} e={e} k={k}: K^0={fmt(k0)}, K^1={fmt(k1)};  dual e'={k} k'={e}: K^0={fmt(d0)}, K^1={fmt(d1)}")
check("K^{k,*}(M_e) = K^{e,*+1}(M_k) for g<=2, 0<=e,k<=5", ok)

# Hopf example of the text: S^3 (e=1,k=0) against S^2 x S^1 (e=0,k=1)
check("S^3 untwisted vs S^2xS^1 twisted by 1", circle_bundle_groups(0, 1, 0) == tuple(reversed(circle_bundle_groups(0, 0, 1))))
check("S^2 x S^1 untwisted: K^1 = Z^2 (not Z)", circle_bundle_groups(0, 0, 0)[1] == (2, []))

# ----------------------------------------------------------------------
# 6. Simply connected closed 6-manifolds with torsion-free cohomology
# ----------------------------------------------------------------------
print("\n6. Twisted K-theory of simply connected six-manifolds (Calabi-Yau threefolds)")


def symplectic_form(b3):
    """Standard unimodular alternating form on Z^{b3} (b3 even)."""
    n = b3 // 2
    M = [[0] * b3 for _ in range(b3)]
    for i in range(n):
        M[i][n + i] = 1
        M[n + i][i] = -1
    return M


def twisted_k_6manifold(b2, b3, omega, Hvec):
    """X simply connected, closed, oriented, H^*(X;Z) torsion-free,
    b_4 = b_2, H^1 = H^5 = 0; omega = intersection form on H^3; Hvec = the
    twist in H^3.  Differentials: d_3 = Sq^3_Z - H cup (-); Sq^3_Z vanishes
    on H^0 and H^2 (degree), and on H^3 it is a Bockstein image, hence
    torsion, hence zero in H^6 = Z.  So d_3(1) = -H and d_3(y) = -omega(H,y)
    on H^3; all other differentials vanish for degree reasons."""
    H = Matrix(Hvec)
    Om = Matrix(omega)
    functional = (H.T * Om)  # row vector y -> omega(H,y)
    # E^0 = ker(Z -> H^3, 1 -> H)
    e0 = (1, []) if all(h == 0 for h in Hvec) else (0, [])
    e2 = (b2, [])
    e4 = (b2, [])
    # E^6 = Z / image of functional
    gcd_f = 0
    for v in functional:
        gcd_f = gcd(gcd_f, abs(int(v)))
    e6 = group_from_relations(1, [[gcd_f]])
    # E^3 = ker(functional) / Z H
    if all(h == 0 for h in Hvec):
        e3 = (b3, [])
    else:
        K = functional.nullspace()           # rational basis of the kernel
        # integral kernel = saturation; compute via Smith form of the functional
        # Kernel of a surjective-up-to-gcd map Z^{b3} -> Z is a free direct summand
        # of rank b3-1.  The quotient by Z H: write H in that lattice.
        # Use the Smith normal form of the 1 x b3 matrix to get a unimodular change
        # of basis U with functional * U = (gcd, 0, ..., 0).
        A = functional
        S, U = _smith_with_transform(A)
        # columns 1..b3-1 of U span the kernel lattice
        Hnew = U.inv() * H                   # coordinates of H in the new basis
        assert Hnew[0] == 0
        rel = [[int(Hnew[i]) for i in range(1, b3)]]
        e3 = group_from_relations(b3 - 1, rel)
    k0 = (e0[0] + e2[0] + e4[0] + e6[0], e6[1])   # all extensions split (free quotients)
    k1 = e3
    return k0, k1, (e0, e2, e3, e4, e6)


def _smith_with_transform(A):
    """Column operations only (A has one row): returns (S, U) with A*U = S,
    S = (g,0,...,0), U unimodular.  Extended Euclid on the row."""
    from sympy import eye
    A = Matrix(A)
    ncol = A.cols
    U = eye(ncol)
    row = list(A)
    # repeatedly reduce
    while True:
        nz = [j for j in range(ncol) if row[j] != 0]
        if len(nz) <= 1:
            break
        j0 = min(nz, key=lambda j: abs(row[j]))
        for j in nz:
            if j == j0:
                continue
            q = row[j] // row[j0]
            row[j] -= q * row[j0]
            U[:, j] = U[:, j] - q * U[:, j0]
    nz = [j for j in range(ncol) if row[j] != 0]
    if nz and nz[0] != 0:
        j = nz[0]
        U.col_swap(0, j)
        row[0], row[j] = row[j], row[0]
    if row[0] < 0:
        U[:, 0] = -U[:, 0]
        row[0] = -row[0]
    S = Matrix([row])
    assert A * U == S
    return S, U


examples = [("quintic", 1, 204), ("bicubic in P^2xP^2", 2, 168), ("Schoen (19,19)", 19, 40)]
for name, b2, b3 in examples:
    om = symplectic_form(b3)
    for m in [0, 1, 2, 3]:
        Hvec = [m] + [0] * (b3 - 1)       # H = m times a primitive class
        k0, k1, pieces = twisted_k_6manifold(b2, b3, om, Hvec)
        print(f"  {name:<30} m={m}: K^0 = {fmt(k0):<14} K^1 = {fmt(k1)}")
        RESULTS.setdefault("six_manifolds", []).append([name, m, fmt(k0), fmt(k1)])
        if m == 0:
            check(f"{name} untwisted ranks 2+2b2 and b3", k0 == (2 + 2 * b2, []) and k1 == (b3, []))
        else:
            check(f"{name} m={m}: K^0 = Z^(2b2)+Z/m, K^1 = Z^(b3-2)+Z/m",
                  k0 == (2 * b2, [m] if m > 1 else []) and k1 == (b3 - 2, [m] if m > 1 else []))
# a non-coordinate primitive class and a non-primitive one on the quintic
om = symplectic_form(204)
Hvec = [3, 5] + [0] * 202      # primitive (gcd 1)
k0, k1, _ = twisted_k_6manifold(1, 204, om, Hvec)
check("quintic, H=(3,5,0,...): K^0 = Z^2, K^1 = Z^202", k0 == (2, []) and k1 == (202, []))
Hvec = [6, 10, 4] + [0] * 201  # = 2 (3,5,2), divisibility 2
k0, k1, _ = twisted_k_6manifold(1, 204, om, Hvec)
check("quintic, H=2(3,5,2,0,...): K^0 = Z^2+Z/2, K^1 = Z^202+Z/2", k0 == (2, [2]) and k1 == (202, [2]))


# ----------------------------------------------------------------------
# 6b. Torsion in H^3: C_H = Tor(H^3 / Z H) and the surgery manifolds X_n
# ----------------------------------------------------------------------
print("\n6b. Six-manifolds with torsion in H^3 (Theorem 1.1 in general, the presentation of C_H and the surgery manifolds X_n)")


def C_H_presentation(m, tors, t):
    """Tor(H^3/ZH) from the presentation of C_H in the paper: generators
    u, c_1..c_r; relations m u + sum t_j c_j = 0 and n_j c_j = 0."""
    r = len(tors)
    rels = [[m] + list(t)] + [[0] + [tors[j] if i == j else 0 for i in range(r)] for j in range(r)]
    return group_from_relations(1 + r, rels)


def twisted_k_6manifold_general(b2, b3, omega, Hfree, tors, Htors):
    """X closed spin^c six-manifold with H_1 = 0; H^3 = Z^b3 (+) T,
    T = (+) Z/n_j (tors), H = (Hfree, Htors).  Returns K^{H,0}, K^{H,1}
    as (rank, torsion) by Theorem 1.1: K^1 = ker f_H / Z H computed
    directly from the lattice, torsion of K^0 = character group of the
    torsion of K^1 (same invariant factors), ranks 2+2b2-2d and b3-2d."""
    r = len(tors)
    tors_rels = lambda off: [[0] * off + [tors[j] if i == j else 0 for i in range(r)] for j in range(r)]
    if all(h == 0 for h in Hfree):             # torsion twist: f_H = 0, K^1 = H^3 / <H>
        rels = [[0] * b3 + list(Htors)] + tors_rels(b3)
        k1 = group_from_relations(b3 + r, rels)
        delta = 0
    else:
        H = Matrix(Hfree)
        functional = (H.T * Matrix(omega))
        S, U = _smith_with_transform(functional)
        Hnew = U.inv() * H
        assert Hnew[0] == 0                      # H lies in ker f_H
        rels = [[int(Hnew[i]) for i in range(1, b3)] + list(Htors)] + tors_rels(b3 - 1)
        k1 = group_from_relations(b3 - 1 + r, rels)
        delta = 1
    k0 = (2 + 2 * b2 - 2 * delta, k1[1])
    return k0, k1


# (a) the cyclic formula Z/g + Z/(mn/g), g = gcd(m,n,t), and the splitting criterion
ok_formula, ok_split, ok_order = True, True, True
for m in range(1, 7):
    for n in range(1, 9):
        for t in range(n):
            g = gcd(gcd(m, n), t)
            C = C_H_presentation(m, [n], [t])
            expected = sorted(d for d in (g, m * n // g) if d > 1)
            ok_formula &= C == (0, expected)
            ok_order &= (lambda rk, tt: rk == 0 and math.prod(tt) == m * n)(*C)
            # splits iff C = Z/n + Z/m up to the invariant factors of that direct sum
            split = group_from_relations(2, [[m, 0], [0, n]])
            ok_split &= (C == split) == (t % gcd(m, n) == 0)
check("C_H = Z/g + Z/(mn/g), g = gcd(m,n,t), for m<=6, n<=8, all t", ok_formula)
check("|C_H| = m n for T = Z/n", ok_order)
check("0 -> T -> C_H -> Z/m -> 0 splits iff gcd(m,n) | t", ok_split)

# (b) the surgery manifolds X_n: b2 = 0, b3 = 2, H^3 = Z a + Z b + Z/n c, H = m a + t c
om2 = [[0, 1], [-1, 0]]
surgery_examples = [((2, 2, 0), [2, 2]), ((2, 2, 1), [4]), ((2, 4, 0), [2, 4]),
                    ((2, 4, 1), [8]), ((3, 6, 3), [3, 6]), ((4, 8, 2), [2, 16])]
for (m, n, t), expected in surgery_examples:
    k0, k1 = twisted_k_6manifold_general(0, 2, om2, [m, 0], [n], [t])
    print(f"  X_{n}, H = {m}a + {t}c: K^0 = {fmt(k0):<12} K^1 = {fmt(k1)}")
    RESULTS.setdefault("surgery_manifolds", []).append([n, m, t, fmt(k0), fmt(k1)])
    check(f"X_{n}, H={m}a+{t}c: both groups finite with invariant factors {expected}",
          k0 == (0, expected) and k1 == (0, expected))
check("X_2: fluxes 2a and 2a+c give Z/2+Z/2 and Z/4",
      twisted_k_6manifold_general(0, 2, om2, [2, 0], [2], [0])[1] == (0, [2, 2]) and
      twisted_k_6manifold_general(0, 2, om2, [2, 0], [2], [1])[1] == (0, [4]))

# (c) the lattice computation agrees with the presentation for non-cyclic T and a
#     flux with several nonzero free coordinates (quintic numbers with torsion added)
om = symplectic_form(8)
agree = True
for tors in ([2, 4], [2, 2, 4], [3, 3], [2, 3, 4]):
    for Hfree, m in (([2, 0, 0, 0, 0, 0, 0, 0], 2), ([6, 10, 4, 0, 0, 0, 0, 0], 2), ([3, 6, 0, 9, 0, 0, 0, 0], 3)):
        for Ht in ([0] * len(tors), [1] * len(tors), [tors[j] - 1 for j in range(len(tors))]):
            k0, k1 = twisted_k_6manifold_general(3, 8, om, Hfree, tors, Ht)
            C = C_H_presentation(m, tors, Ht)
            agree &= k1 == (8 - 2, C[1]) and k0 == (2 + 6 - 2, C[1])
            agree &= math.prod(C[1]) == m * math.prod(tors)
check("lattice computation = presentation of C_H for T in {Z/2+Z/4, Z/2+Z/2+Z/4, Z/3+Z/3, Z/2+Z/3+Z/4}", agree)

# (d) torsion twist and no twist: K^1 = Z^b3 + T/<H>, K^0 = Z^(2+2b2) + (T/<H>)^v
k0, k1 = twisted_k_6manifold_general(3, 8, om, [0] * 8, [2, 4], [0, 0])
check("H = 0 with T = Z/2+Z/4: K^0 = Z^8 + Z/2 + Z/4, K^1 = Z^8 + Z/2 + Z/4", k0 == (8, [2, 4]) and k1 == (8, [2, 4]))
k0, k1 = twisted_k_6manifold_general(3, 8, om, [0] * 8, [2, 4], [0, 2])
check("torsion twist H = (0,2) in Z/2+Z/4: K^1 = Z^8 + Z/2 + Z/2", k1 == (8, [2, 2]) and k0 == (8, [2, 2]))
# the torsion-free function and the general one agree on the quintic
k0a, k1a, _ = twisted_k_6manifold(1, 204, symplectic_form(204), [6, 10, 4] + [0] * 201)
k0b, k1b = twisted_k_6manifold_general(1, 204, symplectic_form(204), [6, 10, 4] + [0] * 201, [], [])
check("general and torsion-free computations agree on the quintic with H = 2(3,5,2,0,...)", (k0a, k1a) == (k0b, k1b))


# ----------------------------------------------------------------------
# 7. Real K-theory torsion charges of Type I (Witten): KO^{-n}(pt)
# ----------------------------------------------------------------------
print("\n7. Type I torsion charges from KO~(S^{9-p}) (Bott periodicity for KO)")
KO = {0: "Z", 1: "Z/2", 2: "Z/2", 3: "0", 4: "Z", 5: "0", 6: "0", 7: "0"}
for p in [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8]:
    n = 9 - p
    print(f"  D{p}: KO~(S^{n}) = KO^{{-{n}}}(pt) = {KO[n % 8]}")
check("Type I D(-1), D0, D7, D8 carry Z/2 charges", all(KO[(9 - p) % 8] == "Z/2" for p in (-1, 0, 7, 8)))
RESULTS["type_I_torsion"] = {f"D{p}": KO[(9 - p) % 8] for p in range(-1, 9)}

print(f"\n{PASSED} checks passed, 0 failed")
RESULTS["checks_passed"] = PASSED
with open("results.json", "w") as f:
    json.dump(RESULTS, f, indent=1)
