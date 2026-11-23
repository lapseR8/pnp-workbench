import Mathlib
import PEqualNP.Defs

theorem PEqualNP.Outside.L_meta_733 : ∀ (n : ℕ) (G : SimpleGraph (Fin n)) (X : Set (Fin n)),
  (∀ u v, G.Adj u v → (u ∈ X ↔ v ∉ X)) →
  ∀ I : Set (Fin n),
    ((∀ u ∈ I, ∀ v ∈ I, ¬ G.Adj u v) ∧ (∀ v, v ∉ I → ∃ u ∈ I, G.Adj v u)) ↔
    ∃ B : Set (Fin n), (∀ b ∈ B, b ∉ X) ∧
      (∀ y, y ∉ X → ((∀ x, G.Adj y x → ∃ b ∈ B, G.Adj x b) ↔ y ∈ B)) ∧
      ∀ v, v ∈ I ↔ (v ∈ B ∨ (v ∈ X ∧ ∀ b ∈ B, ¬ G.Adj v b)) := by
  intro n G X hX I
  constructor
  · rintro ⟨hind, hdom⟩
    refine ⟨{v | v ∈ I ∧ v ∉ X}, fun b hb => hb.2, ?_, ?_⟩
    · intro y hy
      constructor
      · intro h
        by_contra hyB
        have hyI : y ∉ I := fun hyI => hyB ⟨hyI, hy⟩
        obtain ⟨u, huI, hyu⟩ := hdom y hyI
        obtain ⟨b, hb, hub⟩ := h u hyu
        exact hind u huI b hb.1 hub
      · intro hyB x hyx
        exact ⟨y, hyB, hyx.symm⟩
    · intro v
      constructor
      · intro hv
        by_cases hvX : v ∈ X
        · exact Or.inr ⟨hvX, fun b hb => hind v hv b hb.1⟩
        · exact Or.inl ⟨hv, hvX⟩
      · rintro (h | ⟨hvX, h⟩)
        · exact h.1
        · by_contra hvI
          obtain ⟨u, huI, hvu⟩ := hdom v hvI
          have huX : u ∉ X := (hX v u hvu).1 hvX
          exact h u ⟨huI, huX⟩ hvu
  · rintro ⟨B, hB, hsat, hI⟩
    constructor
    · intro u hu v hv huv
      rcases (hI u).1 hu with hu' | ⟨huX, hu'⟩ <;>
        rcases (hI v).1 hv with hv' | ⟨hvX, hv'⟩
      · have := hX u v huv
        exact hB v hv' (by_contra fun h => hB u hu' (this.2 h))
      · exact hv' u hu' huv.symm
      · exact hu' v hv' huv
      · exact (hX u v huv).1 huX hvX
    · intro v hvI
      have hvB : v ∉ B := fun h => hvI ((hI v).2 (Or.inl h))
      by_cases hvX : v ∈ X
      · have : ¬ ∀ b ∈ B, ¬ G.Adj v b := fun h => hvI ((hI v).2 (Or.inr ⟨hvX, h⟩))
        push Not at this
        obtain ⟨b, hb, hvb⟩ := this
        exact ⟨b, (hI b).2 (Or.inl hb), hvb⟩
      · have : ¬ ∀ x, G.Adj v x → ∃ b ∈ B, G.Adj x b := fun h => hvB ((hsat v hvX).1 h)
        push Not at this
        obtain ⟨x, hvx, hx⟩ := this
        have hxX : x ∈ X := by
          by_contra hxX
          exact hvX ((hX v x hvx).2 hxX)
        exact ⟨x, (hI x).2 (Or.inr ⟨hxX, fun b hb => hx b hb⟩), hvx⟩
