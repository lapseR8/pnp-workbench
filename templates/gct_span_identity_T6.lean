import Mathlib

namespace T6

abbrev M := Fin 4 × Fin 5 → ℂ
abbrev aC := M
def aO (M:aC):ℂ:=∑e:Fin 4 ↪ Fin 5,∏r,M (r,e r)
abbrev aHw:=Fin 5×Fin 5×Fin 5×Fin 5
def yH:aHw≃(Fin 4→Fin 5) where
 toFun p:= ![p.1,p.2.1,p.2.2.1,p.2.2.2]
 invFun z:=(z 0,z 1,z 2,z 3)
 left_inv _:=rfl
 right_inv z:=by ext i;fin_cases i<;>rfl
lemma cK (M:aC):
  aO M=∑j,∑l,∑m,∑n,
   if j=l∨j=m∨j=n∨l=m∨l=n∨m=n then 0 else
   M (0,j)*M (1,l)*M (2,m)*M (3,n):=by
 classical
 let u (z:Fin 4→Fin 5):=∏r,M (r,z r)
 have A:aO M=∑z:Fin 4→Fin 5,if Function.Injective z then u z else 0:=by
  let E:=Equiv.subtypeInjectiveEquivEmbedding (Fin 4) (Fin 5)
  change(∑e:Fin 4 ↪ Fin 5,u e)=_;rw[←E.sum_comp (fun e:Fin 4 ↪ Fin 5=>u e),←Finset.sum_filter];exact(Finset.sum_subtype _ (by simp) u).symm
 rw[A,←yH.sum_comp (fun z=>if Function.Injective z then u z else 0)];simp only[Fintype.sum_prod_type];refine sum_congr rfl fun a _=>sum_congr rfl fun b _=>sum_congr rfl fun c _=>sum_congr rfl fun d _=> ?_
 let z:=yH (a,b,c,d)
 have p:z 0=a:=rfl;have q:z 1=b:=rfl;have s:z 2=c:=rfl;have t:z 3=d:=rfl;have f:Function.Injective z↔¬(a=b∨a=c∨a=d∨b=c∨b=d∨c=d):=by
  constructor
  ·intro o;have h (i j:Fin 4) (g:i≠j):z i≠z j:=fun h=>g (o h)
   exact not_or_intro (h 0 1 (by decide)) (not_or_intro (h 0 2 (by decide))
    (not_or_intro (h 0 3 (by decide)) (not_or_intro (h 1 2 (by decide))
      (not_or_intro (h 1 3 (by decide)) (h 2 3 (by decide))))))
  ·intro g i j v;simp only[not_or] at g;fin_cases i<;>fin_cases j<;>first | rfl | grind
 change(if Function.Injective z then u z else 0)=_;by_cases g:a=b∨a=c∨a=d∨b=c∨b=d∨c=d
 ·have k:¬Function.Injective z:=fun h=>f.mp h g
  rw[ite_eq_right k,ite_eq_left g]
 ·rw[ite_eq_left (f.mpr g),ite_eq_right g];change(∏r,M (r,z r))=_;simp[Fin.prod_univ_succ,p,q,s,t,mul_assoc]

def e (i : Fin 4) (j : Fin 5) : M := fun p => if p = (i, j) then 1 else 0

def B : Fin 6 → M := ![e 0 1 - e 0 4, e 0 2 - e 0 4, e 0 3 - e 0 4,
  e 1 0 - e 0 0, e 2 0 - e 0 0, e 3 0 - e 0 0]

set_option maxHeartbeats 4000000 in
theorem ident (c : Fin 6 → ℂ) :
    aO (1 + ∑ k, c k • B k) = 120 := by
  rw [cK]
  simp only [Fin.sum_univ_five, Fin.sum_univ_six, Pi.add_apply, Finset.sum_apply, Pi.smul_apply,
    smul_eq_mul, B, e, Pi.one_apply]
  simp
  ring

theorem lin : LinearIndependent ℂ B := by
  rw [Fintype.linearIndependent_iff]
  intro g hg
  have h := fun p => congrFun hg p
  simp only [Finset.sum_apply, Pi.smul_apply, smul_eq_mul, Fin.sum_univ_six, B, e,
    Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val, Pi.sub_apply,
    Pi.zero_apply] at h
  have h1 := h (0, 1); have h2 := h (0, 2); have h3 := h (0, 3)
  have h4 := h (1, 0); have h5 := h (2, 0); have h6 := h (3, 0)
  simp at h1 h2 h3 h4 h5 h6
  intro i; fin_cases i <;> simp_all

theorem main : ∃ L : Submodule ℂ (Fin 4 × Fin 5 → ℂ),
    (∀ t ∈ L, (∑ f : Fin 4 ↪ Fin 5, ∏ i, (1 + t (i, f i))) = 120) ∧ Module.finrank ℂ L = 6 := by
  refine ⟨Submodule.span ℂ (Set.range B), ?_, ?_⟩
  · intro t ht
    obtain ⟨c, rfl⟩ := (Submodule.mem_span_range_iff_exists_fun ℂ).1 ht
    exact ident c
  · rw [finrank_span_eq_card lin]; simp

end T6
