import Mathlib
import PEqualNP.Defs
import PEqualNP.MetaComplexity.L104
import PEqualNP.MetaComplexity.L105

namespace PEqualNP.Outside.L_meta_748

def starG : SimpleGraph (Fin 6) where
  Adj i j := (i = 0 ∧ j ≠ 0) ∨ (j = 0 ∧ i ≠ 0)
  symm := ⟨fun i j h => by rcases h with h | h <;> [exact Or.inr h; exact Or.inl h]⟩
  loopless := ⟨fun i h => by rcases h with h | h <;> exact h.2 h.1⟩

def starDec : DecidableRel starG.Adj := fun i j =>
  (inferInstance : Decidable ((i = 0 ∧ j ≠ 0) ∨ (j = 0 ∧ i ≠ 0)))

def circ : PEqualNP.CircuitLowerBounds.Circuit :=
  ⟨[⟨PEqualNP.CircuitLowerBounds.GateOp.or, [1, 2]⟩,
    ⟨PEqualNP.CircuitLowerBounds.GateOp.or, [6, 3]⟩,
    ⟨PEqualNP.CircuitLowerBounds.GateOp.or, [7, 4]⟩,
    ⟨PEqualNP.CircuitLowerBounds.GateOp.or, [8, 5]⟩,
    ⟨PEqualNP.CircuitLowerBounds.GateOp.and, [0, 9]⟩], 10⟩

theorem circ_eval : ∀ x : Fin 6 → Bool, circ.eval (List.ofFn x) =
    ((x 0 && x 1) || (x 0 && x 2) || (x 0 && x 3) || (x 0 && x 4) || (x 0 && x 5)) := by
  decide

theorem upper (s : ℕ) (hs : 5 ≤ s) : PEqualNP.MetaComplexity.HasCircuitSize 6
    (fun x : Fin 6 → Bool => (x 0 && x 1) || (x 0 && x 2) || (x 0 && x 3) || (x 0 && x 4) ||
      (x 0 && x 5)) s := by
  refine ⟨circ, ?_, ?_, circ_eval⟩
  · intro g hg
    simp only [circ, List.mem_cons, List.not_mem_nil, or_false] at hg
    rcases hg with rfl | rfl | rfl | rfl | rfl <;> simp
  · show 5 ≤ s
    exact hs

theorem connected : starG.Connected := by
  have h0 : ∀ v : Fin 6, starG.Reachable 0 v := by
    intro v
    by_cases hv : v = 0
    · subst hv; exact SimpleGraph.Reachable.refl _
    · exact SimpleGraph.Adj.reachable (Or.inl ⟨rfl, hv⟩)
  exact ⟨fun u v => (h0 u).symm.trans (h0 v)⟩

theorem acyclic : starG.IsAcyclic := by
  let _ := starDec
  have hcard : Nat.card starG.edgeSet + 1 = Nat.card (Fin 6) := by
    rw [Nat.card_eq_fintype_card, Nat.card_eq_fintype_card, ← SimpleGraph.edgeFinset_card]
    decide
  exact ((SimpleGraph.isTree_iff_connected_and_card).2 ⟨connected, hcard⟩).isAcyclic

theorem fun_eq : ∀ x : Fin 6 → Bool,
    ((x 0 && x 1) || (x 0 && x 2) || (x 0 && x 3) || (x 0 && x 4) || (x 0 && x 5)) =
      @decide (∃ i j : Fin 6, starG.Adj i j ∧ x i = true ∧ x j = true)
        (letI := starDec; inferInstance) := by
  decide

theorem lower (s : ℕ) (h : PEqualNP.MetaComplexity.HasCircuitSize 6
    (fun x : Fin 6 → Bool => (x 0 && x 1) || (x 0 && x 2) || (x 0 && x 3) || (x 0 && x 4) ||
      (x 0 && x 5)) s) : 5 ≤ s := by
  obtain ⟨C, hC, hs, hf⟩ := h
  let _ := starDec
  have key := PEqualNP.MetaComplexity.L104 6 1 s starG acyclic
    ⟨0, 1, Or.inl ⟨rfl, by decide⟩⟩
    (fun T hT => by
      rcases hT 0 1 (Or.inl ⟨rfl, by decide⟩) with h | h <;> exact Finset.card_pos.2 ⟨_, h⟩)
    ⟨C, hC, hs, fun x => (hf x).trans (fun_eq x)⟩
  have he : starG.edgeFinset.card = 5 := by decide
  omega

end PEqualNP.Outside.L_meta_748

theorem PEqualNP.Outside.L_meta_748 : ∀ s : ℕ, PEqualNP.MetaComplexity.HasCircuitSize 6 (fun x : Fin 6 → Bool => (x 0 && x 1) || (x 0 && x 2) || (x 0 && x 3) || (x 0 && x 4) || (x 0 && x 5)) s ↔ 5 ≤ s := by
  intro s
  exact ⟨PEqualNP.Outside.L_meta_748.lower s, PEqualNP.Outside.L_meta_748.upper s⟩
