import Mathlib
import PEqualNP.Defs

open PEqualNP.CircuitLowerBounds in
def C8 : Circuit := ⟨[⟨.or,[2,4]⟩,⟨.or,[2,5]⟩,⟨.and,[1,9]⟩,⟨.or,[1,7]⟩,⟨.or,[4,11]⟩,⟨.or,[5,7]⟩,⟨.and,[6,13]⟩,⟨.and,[0,12]⟩,⟨.and,[3,8]⟩,⟨.or,[10,16]⟩,⟨.or,[14,17]⟩,⟨.or,[15,18]⟩], 19⟩

set_option maxRecDepth 100000 in
set_option maxHeartbeats 4000000 in
theorem t8 : ∀ x : Fin 8 → Bool, (C8.eval (List.ofFn x) = true ↔ ∃ u v : Fin 8, (u.val, v.val) ∈ ([(0, 1), (1, 2), (2, 3), (3, 4), (4, 0), (1, 5), (5, 6), (6, 7), (7, 0)] : List (ℕ × ℕ)) ∧ x u = true ∧ x v = true) := by
  decide
