import Std

/-!
Dependency-free structural core for FIN Programs P379--P393.

The separate generated arithmetic-reflection file checks the large terminal
integer predicates.  This file records only structural implications.
-/

namespace FINPrograms379To393

theorem damping_completion_square
    (legacy strict completion : Nat → Rat)
    (identity :
      ∀ d, completion d * legacy d = strict d) :
    ∀ d, completion d * legacy d = strict d :=
  identity

theorem pullback_radial_naturality
    {Vertex Distance Value : Type}
    (distanceSource distanceTarget : Vertex → Vertex → Distance)
    (completion : Distance → Value)
    (map : Vertex → Vertex)
    (isometric :
      ∀ x y,
        distanceTarget (map x) (map y) = distanceSource x y) :
    ∀ x y,
      completion (distanceTarget (map x) (map y)) =
      completion (distanceSource x y) := by
  intro x y
  rw [isometric x y]

abbrev Vec5 := Fin 5 → Nat

def generator (selected : Fin 5) : Vec5 :=
  fun index => if index = selected then 1 else 0

theorem generator_has_no_other_coordinate
    (source target : Fin 5)
    (different : source ≠ target) :
    generator source target = 0 := by
  simp [generator, Ne.symm different]

theorem additive_cancellation
    (a b c : Nat)
    (bound : a + c ≥ b + c) :
    a ≥ b := by
  exact Nat.le_of_add_le_add_right bound

theorem positive_discrete_scale_not_half_invariant
    (rho : Nat)
    (positive : rho > 0) :
    rho ≠ rho / 2 := by
  omega

theorem realization_nonuniqueness_from_same_grade
    {Grade Realization : Type}
    (grade : Realization → Grade)
    (left right : Realization)
    (same : grade left = grade right)
    (different : left ≠ right) :
    ∃ x y, grade x = grade y ∧ x ≠ y :=
  ⟨left, right, same, different⟩

end FINPrograms379To393
