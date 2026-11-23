# L-gct-179 algebra notes (attacker C, exact algebra)  [code: scratchpad/alg/]

## C1. E2 exactly (sympy)
Hessian of E2 on C^20: eigenvalues 12 (J), -4 (x3, row-type u(x)1), -3 (x4, col-type 1(x)v), 1 (x12, doubly centred).
rank 20 on C^20, rank 19 on H={E1=0}. So Q nondegenerate on H: complex max isotropic = 9, real signature (12,7) -> real bound 7.
Constraint set used everywhere: for every multiset {a1..ak} of basis indices, full polarization P_k(b_a1..b_ak)=0, k=1..4.

## C2. NEW 6-dim family (not row∪col!) found by LM seeded near the cross, then fitted exactly (sympy verified, all polarizations 0)
F_b (support rows{0,2} ∪ cols{0,2}), parameter b, a=b/(1-b), c=3b-1:
  v1=e00-e10, v2=e30-e10, v3=col2 = (-1/3,-1/3,1,-1/3)^T,
  v4=row0 (0,1,0,-a,a-1), v5=row2 (1,-b,0,b-1,0), v6=row2 (0,-c,0,c-1,1).
  Conditions on (a,b,c) from E1..E4 are exactly {ab-a+b=0, ac-2a+c+1=0}.  So V6 has positive-dim components
  (rigidity-of-6-spaces route is dead).  Cross spaces themselves are numerically isolated (LM from eps=0.05 perturbations returns to them).
## C3. Maximality, exact: for cross (both variants a: row-sum-off-c=0, b: col-sum-off-r=0) and for F_b (b=1/5,2/7,-3/11):
  V1 := {v : all polarizations linear in v vanish on L6} = L6 ⊕ <w>, w ≡ J - 20 e_rc mod L6 (rc = cross centre; (2,2) for F_b),
  and E2(w)=-120, E3(w)=-480, E4(w)=-360, all mixed terms 0 -> Groebner basis (c^2). None of these extends to 7.
  Mechanism: dE_k(x)[J] ∝ E_{k-1}(x)=0 (k>=2) and dE_k(x)[e_rc] = (k-1)-minor sums of x with row r, col c deleted = 0 on L6.
  Zariski tangent space of the Fano scheme V6 at cross: dim 27 (var a), 28 (var b) mod p — very non-reduced, useless for rigidity.

## C4. F_b is maximal for EVERY parameter value, incl. the closure (exact; scripts alg/fbsym.py, alg/fbsym2.py)
- Solution set {ab-a+b=0, ac-2a+c+1=0} is ONE irreducible rational curve (b=1 and a=-1 are inconsistent): a=b/(1-b), c=3b-1.
  Homogeneous basis polynomial in b (no poles): v1=e00-e10, v2=e30-e10, v3=col2 (-1,-1,3,-1)^T,
  v4=row0 (0,1-b,0,-b,2b-1), v5=row2 (1,-b,0,b-1,0), v6=row2 (-3,1,0,1,1). Independent for every b in C and at b=inf
  (limit: v4=(0,-1,0,-1,2), v5=(0,-1,0,1,0)); sympy: all E_k vanish identically in b, and at b=inf. So the family is a P^1.
- Linear system for V1 (all polarizations linear in v): matrix A(b) (16x20, entries deg<=2 in b). gcd over Z[b] of
  det(R A(b) C) for 5 random integer R (13xN), C (20x13) is the constant 12  =>  rank A(b)=13 for EVERY complex b; rank 13 at b=inf.
  Hence dim V1 = 7 on all of P^1.  Explicitly V1 = F_b ⊕ <w>, w = J - 20 e_22 (for all b): dE_k(x)[J] ∝ E_{k-1}(x)=0, and
  dE_k(x)[e_22] = (k-1)-minor sums of x with row2,col2 deleted, which is a zero-sum row0∪col0 cross pattern -> 0.
  For v = x + s w (x in F_b): E2(v) = s^2 E2(w) = -120 s^2 (beta(x,w)=0). So NO extension anywhere on the closure of F_b.
- F_b vs A's families: every member of F_b ∪ {b=inf} has column-pure dims (col0:2, col2:1, others 0) — checked by hand for all b —
  so it never lies in the S4xS5-orbit of W6 (which has a 3-dim zero-sum column-pure part); and V_ζ = C^4⊗φ5 ⊄ F_b since rows 1,3 of
  F_b vanish on cols 1,3,4. The closure of F_b does not meet these families.
- Second numeric 6-dim solution (LM with support rows{0,2}∪cols{0,2}; u_2=2, column-pure col0:1, col2:1; not F_b, not a cross):
  numerically dim V1 = 7 too, extra direction has E2 ≈ -1.92-0.57i ≠ 0. Same "maximal at first order" pattern.

## C5. Equivariant exact search (L invariant under a prime-order g in S4xS5)  [alg/equiv.py, alg/gen_sing.py]
Eigen-decomposition C^20 = ⊕ V_χ (orbit Fourier basis over F_30091, q≡1 mod 15); L = ⊕ L_χ, constraints d_χ+d_{-χ} <= dim V_χ.
dim V: c5 (col 5-cycle) (4,4,4,4,4); c3 (col 3-cycle) (12,4,4); r3 (row 3-cycle) (10,5,5); r3c3 (8,6,6).
Local sympy run, COORDINATE charts (complete for those dim-vectors), c5, d=7: dim-vectors (d_0..d_4)=(0,0,0,3,4),(0,0,0,4,3),(0,0,1,2,4),
(0,0,1,3,3),(0,0,2,1,4),(0,0,2,2,3),(0,0,3,1,3): ALL charts give unit ideal (464 charts).  Rest moved to Singular on GitHub
(lapseR8/pnp-workbench, workflow gct179alg, run 45266023749): all dim-vectors for c5,c3,r3,r3c3 at d=7 with 2 random generic charts each
(probabilistic genericity), plus sanity runs d=6 (r3c3 must be NONUNIT at (2,2,2): cross at (3,4) is r3c3-invariant).

## C6. "Centre" mechanism — common to EVERY 6-dim space found (crosses, F_b, and the 2nd numeric family: J-20e_22 ∈ V1, resid 3e-12)
Def: cell (r,c) is a centre of a valid L if for all x in L the 3x4 block x' (delete row r, col c) has E1(x')=E2(x')=E3(x')=0.
Lemma (exact, elementary): centre (r,c) ⇒ w := J - 20 e_rc ∈ V1(L), and E2(w)=-120, E3(w)=-480, E4(w)=-360.
  (dE_k(x)[J] = c_k E_{k-1}(x) = 0 for k>=2; dE_k(x)[e_rc] = E_{k-1}(x'); E1(w)=0.)
Corollary: if L7 is a valid 7-dim space, then EVERY hyperplane L' ⊂ L7 (a valid 6-space with V1(L') ⊇ L7) has
  dim V1(L') >= 8 or has NO centre.  [If dim V1(L')=7 and w ∈ V1(L') then V1(L')=L7 ∋ w, but E2(w)=-120≠0.]
So a 7-dim witness needs a P^6 of valid 6-spaces that are all either centre-free or first-order non-maximal; none of the
6-spaces found so far (crosses, F_b on all of P^1, numeric 2nd family) is of either kind.
C6'. Unification with A's §7 profile: the obstructing vector is always w_S = J - λ_S 1_S with S a cell (λ=20), a column (λ=5) or a
  row (λ=4); E1(w_S)=0 and per(J+μ w_S) = 120(1+μ)^3(1-3μ) for all three (direct count). Lemma holds verbatim with
  "centre S": dE_k(x)[1_S]=0 for x∈L, k>=2. Crosses/F_b/2nd numeric family: S = cell; W6: S = column; V_ζ: S = row.
  Corollary of C6 generalises: every hyperplane of a 7-dim witness must be centre-free (all three kinds) or have dim V1 >= 8.
