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
