# L-gct-179 structure notes (attacker A) — live file

Notation: rows 0..3, cols 0..4. u_i = dim(L ∩ row-i-supported), r_i = dim(row_i(L)), p_ij = dim proj to rows i,j.

## 0. Basic facts (from 215 file, valid for ANY L with per(J+t)=120 on L, no dim assumption)
- lO: for every row i, there is a functional f_i on C^5 with f_i(1)=1 and f_i(row_i t)=0 ∀t∈L. So r_i ≤ 4.
- lB: p_ij ≤ 7 (dim of proj to two rows).
- aG7 (GENERAL): if u_0 ≥ 4 then dim L ≤ 7. Mechanism: u_0≥4 forces row_0(L) = {row-sum 0} (dim 4) and
  L = R ⊕ P with R = all row-0 zero-sum vectors, P = L ∩ {row 0 = 0}; condition becomes: for t∈P,
  every 3x4 permanent of rows 1..3 of (J+t) with column j deleted equals 24 (∀j). Then xO gives dim P ≤ 3.

## 1. Ansatz U4: u_0 = 4, dim P = 3  (status: being checked, see §2)
P = 3-dim space of 3x5 matrices (rows x,y,z) with
  - column sums zero: x_l+y_l+z_l=0 (forced)
  - write x=a+b, y=ωa+ω²b, z=ω²a+ωb  (ω primitive cube root of 1)
  - quadratic: Σ a_l b_l - A B + a_j B + A b_j - 2 a_j b_j = 0 ∀j  (A=Σa, B=Σb)
  - cubic: e3(a without coord j) + e3(b without coord j) = 0 ∀j
xO case split: A,B both ≢0 on P ⇒ dim P ≤ 2. So WLOG B≡0 on P.
  - b≡0: then each coordinate a_j is a root of x^3-e1x^2+e2x-e3 (of a itself); a 3-dim linear space forces
    a = 3 equal-value blocks, which needs block sizes 1,1,1 (impossible with 5 coords) ⇒ dim ≤ 2.
  - S = supp(b) proper, |S|=s∈{2,3,4}: a_i = A/2 on S, b_i=0 off S, Σb=0. Quadratic auto. Remaining: cubic.
    Ambient W_S is 4-dim; need a hyperplane of W_S killing 5 cubics. COMPUTED (sympy, all 4 charts, s=2,3,4): NO solution.
  - s=1 forces b=0; s=5 forces a≡0 and reduces to the b≡0 argument.
  - Translation (a,b) <-> per-conditions verified symbolically (coefficients of the 3x4 perms = 6*quad, 6*cubic).
  => **THEOREM (A, new): u_i ≥ 4 ⇒ dim L ≤ 6.** Hence any 7-dim L has u_i ≤ 3 for all rows i. Ansatz U4 is DEAD.

## Orchestrator notes (hand derivations, 20:05 UTC; please double-check)
1. Row+column 6-dim family F(r,c) is MAXIMAL. With r=c=0, L6 = {row0 off col0 sum 0, t00 = -(col0 off row0)}.
   First-order terms of E2, E3 vanishing on all of L6 force v (mod L6) = s*(ones on rows0..3 x cols1..4) + w*e00 with v_i0=0 (i>=1);
   E1(v)=0 gives w=-16s, and then E2(v) = 72 - 192 = -120 != 0. So no 7-dim space contains any F(r,c).
2. Supports covered by 2 lines (2 rows: <=5; 2 columns: <=4; row+column: 6) kill E3/E4 structurally but give at most 6.
   A 7-dim witness needs genuine cancellation in E3 and/or E4 (complex, cf. L-gct-158 real case impossible).
3. On H (E1=0), with t = r(x)1 + 1(x)c + d (r, c sum-zero, d doubly centred): 2*E2 = ||d||^2 - 20||r||^2 - 12||c||^2 (bilinear, no conjugation). Nondegenerate on 19 dims, so E2 alone gives dim <= 9.

## 2. General facts for a 7-dim L (no hfix) extracted from the 215 file  [attacker A, ~T+50min]
Valid for every L with per(J+t)=120 on L and dim L = 7:
- (cT, gM) E1=0 and L is totally isotropic for B (polar of E2); B nondegenerate on the 19-dim {E1=0}.
- (aEq) E3-polar: kK(lP T, X)=0 for T,X∈L  (lP T is the quadratic "gradient" of E3).
- (lO) r_i ≤ 4 (a functional f_i with f_i(1)=1 kills row_i(L)).  (lB) p_ij ≤ 7.
- (aBu, NO hfix needed) 7 + r_i ≤ p_ij + p_ik + u_l for {i,j,k,l} = {0,1,2,3}.
- trivial: Σ_i u_i ≤ 7 (row-supported parts are independent); u_i ≤ r_i.
- (A, new, §1) u_i ≤ 3 for every row.
- (eI = aEw/aEv/aDh/vU/std/mK chain) NO 7-dim L with two identical COLUMNS.
- (uZ, line ~6115) NO 7-dim L with two identical ROWS (row0 t = row1 t ∀t).
- Exact row-pure relation (generalises aG7): for r in U_i := row_i(aZ L i) ⊂ 1^⊥ and ANY t∈L,
  Σ_j r_j g_j(t) = 0, where g_j(t) = per_{3x4} of (J+t) with row i, column j deleted.  I.e. g(t) ∈ U_i^⊥ ∀t∈L.
- Exact column-pure relation: for c in C_j := col_j(L ∩ col-j-supported) ⊂ 1^⊥ (ℂ^4) and any t∈L,
  Σ_i c_i h_i(t) = 0, h_i(t) = per_{3x4} of (J+t) minus row i, col j.
- Different rows' pure parts: U_i ⊥ U_k (bilinear dot) for i≠k (from E2); triple/quadruple conditions from E3/E4.
- The dim-8 lemmas MidCase/mI/HighPair/vK all use dim=8 essentially (they produce a 7-dim subspace with identical
  rows/columns, then kill it with eI/uZ). They do NOT give dim-7 constraints beyond the list above.

## 3. Rigidity of the 6-dim example (confirms orchestrator note 1, computed with sympy)
Tangent-extension space T(L6) = {v : ∇E_k(t)·v = 0 ∀t∈L6, k=1..4} is exactly 7-dim = L6 ⊕ <v*>,
v* = (-16,4,0,0,0 / 0,1,1,1,1 / 0,1,1,1,1 / 0,1,1,1,1). On t+s v*: E1=0, E2=-120 s², E3=-480 s³, E4=-360 s⁴
(all mixed terms vanish; only E_k(v*)≠0). So F(r,c) is maximal.

## 4. Small-shape experiments (LM least squares, random complex starts) — CAUTION
- 3x3: d=3 found (2/3), d=4 never (0/8).   3x4: d=4 found (2/3), d=5 never (0/12, resid floor ~0.09).
- 3x5: d=5 found (2/3), d=6 never (0/8).    4x4: d=6 never (0/6).
- Pattern: max = m+n-3 (= row∪col count) for every tested shape. Suggests 4x5 max = 6 (statement false).
- BUT my solver found 0/3 for 4x5 d=6 (which EXISTS), so random-start LM is unreliable at 4x5; the
  negative results for 4x5 d=7 (0/6, resid ≥0.25) are NOT evidence. Numeric attacker: seed with structure.
- Found 3x4 d=4 solutions that are NOT row∪col type, e.g. one with row 0 of rank 1 (row_0(t)=λ(t)·w0 for a
  fixed w0, sum w0 ≠ 0) and rows 1,2 full support → genuine E3 cancellation exists in max-dim solutions.
  (Solutions are singular/slowly convergent, accuracy only ~1e-3 on entries.)

## 5. Remaining ansatz families for a 7-dim witness (attacker A, ~T+80min)
Since u_i ≤ 3, the natural "largest pure part" ansatz is u_0 = 3:
**Ansatz U3.** L ⊇ e_0 ⊗ U with U = 1^⊥ ∩ φ^⊥ (φ ∈ C^5 a parameter, 3-dim U). Then K := π_{rows1..3}(L) is 4-dim and
- row_0(L) = ker f (f(1)=1) ⊇ U; generically U = ker f ∩ 1^⊥ and row 0 of t = λ(k)ψ (mod U), ψ ∈ ker f \ U, λ ∈ K*.
- Exact equations (k ∈ K, scaling k→ck, S1/Q2/G3 as in the 215 file, applied to the 3 rows of k, vectors in C^5):
  (i) S1(k), Q2(k), G3(k) ∈ span(1, φ)      [g(k) ⊥ U]
  (ii) c²: 2ΣQ2 + 6λ ψ·S1 = 0;  c³: ΣG3 + 2λ ψ·Q2 = 0;  c⁴: λ ψ·G3 = 0  (c¹ is automatic from E1).
  For the row∪col example: φ = e_c, K = column-c part (3-dim, λψ carries the corner). Need K 4-dim for dim 7.
- I am running a numeric LM on U3 (φ free) for d=6 (validation) and d=7; results below when done.
**Ansatz R1 (seen in 3x4 max-dim solutions):** some row has rank 1: row_0(t) = λ(t)·w0 (fixed w0, Σw0 ≠ 0).
Then K = ker λ (dim 6) lives on rows 1..3 with Σ_j g_j(k) = 120 only (weak), and the λ-direction couples.
**Necessary for any witness:** support of a generic element must contain 3- and 4-matchings with genuine
cancellation in E3/E4 (supports covered by 2 lines give ≤ 6); no two identical rows/columns; u_i ≤ 3.
Symmetry suggestions (heuristic only — a finite group need not fix a point of the witness variety):
 stabiliser of a column c and a row r (S3×S4), and Z/5 cyclic on columns (eigenspaces e^{2πi jk/5}) — the U3
 ansatz with φ = e_c is S4-symmetric on the other columns; with φ = (1,ζ,ζ²,ζ³,ζ⁴) it is Z/5-covariant.

**Ansatz V3 (column analogue of U4, NOT covered by the 215 file).** v_c := dim(L ∩ col-c-supported) ≤ 3; if v_c = 3
then L ⊇ all zero-sum col-c vectors, K := π_{cols≠c}(L) is a 4-dim space of 4x4 matrices, and EXACTLY:
  h_i(k) := per_{3x4}((J+k) minus row i) is independent of i (=: s(k)),  per_4(J+k) + s(k)(4+σ(k)) = 120,
  σ(k) = column-c sum = -total(k).  Degree split: row sums of k all equal (ρ, σ=-4ρ); E2,E3 of row-deleted
  3x4 blocks independent of i;  2E2+8E2'+6E1'σ=0, E3+4E3'+2E2'σ=0, E4+E3'σ=0.
  On K0={ρ=0} (dim ≥3, rows x_a ∈ 1^⊥): Gram g_ab=x_a·x_b must satisfy Σ_{b≠a} g_ab = 0 ∀a (bilinearly on K0).
  Row∪col example = V3 with K 3-dim. A witness of type V3 needs K 4-dim. (numeric test running)

## 6. More results (attacker A, ~T+115min) — new families, all RIGID
**(a) New 4-dim solution V_ζ = C^4 ⊗ ζ**, ζ=(1,ζ5,ζ5²,ζ5³,ζ5⁴): E_k(x ζ^T) = k!·e_k(x)·e_k(ζ) = 0 for k≤4
  (e_1..e_4 of the 5th roots of unity vanish). Tangent-extension T(V_ζ) is 7-dim = V_ζ ⊕ {r⊗1 : Σr=0}
  (row-constant part), which is NOT B-isotropic. Hence **no 7-dim L contains C^4⊗ζ** (any L ⊇ V must lie in T(V)),
  and in fact V_ζ extends to at most dim 5 (isotropic part of {r⊗1} is 1-dim).
**(b) New 6-dim family W6 (found by the V3-ansatz LM, then identified exactly):**
  W6 = {zero-sum vectors in column c} ⊕ {x ⊗ w : x supported on 3 rows}, w on the other 4 columns with
  p1=p2=p3=0, i.e. w = scaled 4th roots of unity, e.g. w=(0,-1,1,-i,i). Here E_k(x w^T)=k! e_k(x) e_k(w),
  e_1=e_2=e_3=0, e_4(w)=-1, so need e_4(x)=0 ⇒ x on 3 rows; cross terms with column c vanish.
  NOT of row∪col type (genuine 3-/4-matchings, cancellation via roots of unity). Invariants: v_c=3, u=1 on 3 rows.
  T(W6) = W6 ⊕ <v>, v = (col c ≡ 1, other cols ≡ -1/4) — the non-isotropic "column type" vector. **W6 is rigid.**
**(c) Numerics with the V3 ansatz (L ⊇ zero-sum column c, analytic-Jacobian LM, good conditioning):**
  d=6: 4/8 converge quadratically to 1e-13 (all of W6 type); d=7: 0/10 (resid ≥ 0.25). Good evidence: **v_c = 3 ⇒ dim ≤ 6.**
  U3 ansatz (finite-diff LM) gave nothing at d=6 or 7 — inconclusive (solver weak).
**Pattern:** every maximal family found (row∪col F(r,c), W6) is rigid: its first-order extension space adds exactly
  one vector, and that vector is a pure non-isotropic type (row-constant / column-constant / J-block), failing E2.

## VERDICT (attacker A): unresolved, leaning strongly to "statement_seems_false" (max = 6).

## 7. F_b (attacker C's family) checked by attacker A  [20:44 UTC]
- Verified SYMBOLICALLY in b (sympy, generic b): E1..E4 vanish identically on F_b.
- Invariants (numeric ranks, generic b and b=1/5, 2/7, -3/11): u = (1,0,2,0), v = (2,0,1,0,0),
  r = (3,2,3,2), p01=4,p02=5,p03=4,p12=4,p13=3,p23=4. Support rows{0,2} ∪ cols{0,2} (4 lines).
- So F_b is a NEW TYPE: not U4 (killed), not V3 (max v_c = 2), not U3 (max u_i = 2). It does not contradict
  "v_c = 3 ⇒ dim ≤ 6"; it simply lies outside the ansätze I tested.
- MAXIMAL FOR ALL b (symbolic): T(F_b) = F_b ⊕ <w_b>,
  w_b = [[3,(3b-2)/(2b-1),0,(3b-1)/(2b-1),0],[0,1,0,1,1],[4,0,-16,0,0],[0,1,0,1,1]].
  On F_b ⊕ <w_b>, ALL mixed terms vanish and E2 = -120 s², E3 = -480 s³, E4 = -360 s⁴, independent of b.
  The special points b = 1/2, 1/3, 0, 2/3 were checked numerically: still dim T = 7 and still obstructed.
- **Universal obstruction (observed in all four maximal families: L6, W6, F_b, and the V_ζ extension):**
  T(L)/L is 1-dimensional, and the extra direction w (canonical mod L, because mixed terms vanish) always has
  per(J + μw) = 120 (1+μ)³ (1-3μ) after rescaling μ, i.e. E3²/(E2·E4) = 16/3 and E2²/E4 = -40.
  Example: w = (column c ≡ 1, other columns ≡ -1/4). No witness can contain ANY of these 6-dim spaces.

## 8. HONEST CASE FOR "NO 7-DIM SPACE" (for the release note)  [attacker A, 20:46 UTC]
PROVED (Lean or exact computation):
 P1. dim ≤ 7 (L-gct-215, Lean). On {E1=0}, E2 is a nondegenerate quadratic form, so E2 alone gives dim ≤ 9.
 P2. 215 lemmas valid for every 7-dim L (no hfix): r_i ≤ 4; p_ij ≤ 7; 7 + r_i ≤ p_ij + p_ik + u_l;
     L is B-isotropic and satisfies the E3-polar identity; no two identical columns (eI); no two identical rows (uZ).
 P3. (A, exact) u_i ≥ 4 ⇒ dim ≤ 6. Proof: aG7's reduction, the xO case split, and an exhaustive sympy hyperplane
     search (§1). So every 7-dim L has all row-pure parts of dim ≤ 3.
 P4. (exact tangent computations) Each known maximal-dimension family is non-extendable: F(r,c) (row∪col),
     W6 (column ⊕ roots-of-unity rank-1 block), F_b (all b), and C^4⊗ζ (it cannot even reach dim 6).
     Any 7-dim L containing a 6-dim solution L' must lie in T(L'), and in every case T(L') fails E2–E4.
NUMERICAL EVIDENCE ONLY:
 N1. Smaller shapes (3x3, 3x4, 3x5, 4x4): the maximum found is m+n-3 every time, and one more dimension was never
     found. 4x5 analogue: m+n-3 = 6.
 N2. Attacker B: d=7 objective floors at ≈6e-4 over ~2500 restarts, while d=6 converges to 1e-13.
 N3. V3 ansatz (L ⊇ zero-sum column): d=6 found 4/8 with quadratic convergence; d=7 0/10 (floor ≥ 0.25).
 N4. Universal obstruction profile per(J+μw) = 120(1+μ)³(1-3μ) in all four rigid families.
NOT PROVED: there is no case-complete argument for u_i ≤ 3 or for 6-dim spaces in general; the 8→7 machinery of
 the 215 file does not transfer (its dim-8 lemmas use dim = 8 essentially). Equivalent open statement:
 "for every 6-dim solution L', T(L') is not a solution".
No route to dim 7 seen.
