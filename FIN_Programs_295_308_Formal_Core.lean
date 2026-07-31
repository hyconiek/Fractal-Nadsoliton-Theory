import Std

/-!
FIN Programs P295--P308: dependency-free formal core.

This file machine-checks two general structural results without importing
physical interpretation.

P296:
  staged minimization is a valid joint minimization whenever the inner and
  outer minimizers satisfy their declared universal inequalities.  Positive
  quadratic Schur elimination is an instance of this principle.

P308:
  a free nontrivial action has no globally invariant point, while two
  equivariant maps on a transitive torsor are equal once they agree at one
  supplied point.

The file does not derive the minimizers, a physical point, a selector, a unit,
or a legacy-to-strict bridge.
-/

namespace FINPrograms295To308

theorem staged_minimizer_is_joint
    {X Y Z V : Type}
    (le : V → V → Prop)
    (leTrans : ∀ a b c, le a b → le b c → le a c)
    (F : X → Y → Z → V)
    (zStar : X → Y → Z)
    (yStar : X → Y)
    (innerMinimal :
      ∀ x y z, le (F x y (zStar x y)) (F x y z))
    (outerMinimal :
      ∀ x y, le (F x (yStar x) (zStar x (yStar x)))
        (F x y (zStar x y))) :
    ∀ x y z,
      le (F x (yStar x) (zStar x (yStar x))) (F x y z) := by
  intro x y z
  exact leTrans _ _ _ (outerMinimal x y) (innerMinimal x y z)

theorem attained_minimum_value_unique
    {X Y Z V : Type}
    (le : V → V → Prop)
    (leAntisymm : ∀ a b, le a b → le b a → a = b)
    (F : X → Y → Z → V)
    (zStar : X → Y → Z)
    (yStar : X → Y)
    (x : X)
    (v : V)
    (vLower : ∀ y z, le v (F x y z))
    (vAttained :
      le (F x (yStar x) (zStar x (yStar x))) v) :
    F x (yStar x) (zStar x (yStar x)) = v := by
  exact leAntisymm _ _ vAttained
    (vLower (yStar x) (zStar x (yStar x)))

theorem no_invariant_point_for_free_nontrivial_action
    {G T : Type}
    [One G]
    (act : G → T → T)
    (freeAction : ∀ g t, act g t = t → g = 1)
    (nontrivial : ∃ g : G, g ≠ 1) :
    ¬ ∃ t : T, ∀ g : G, act g t = t := by
  intro invariantPoint
  obtain ⟨t, fixed⟩ := invariantPoint
  obtain ⟨g, hg⟩ := nontrivial
  exact hg (freeAction g t (fixed g))

theorem equivariant_maps_equal_when_pointed
    {G T U : Type}
    (actT : G → T → T)
    (actU : G → U → U)
    (t0 : T)
    (transitiveFromPoint : ∀ t : T, ∃ g : G, actT g t0 = t)
    (f h : T → U)
    (fEquivariant : ∀ g t, f (actT g t) = actU g (f t))
    (hEquivariant : ∀ g t, h (actT g t) = actU g (h t))
    (agreeAtPoint : f t0 = h t0) :
    f = h := by
  funext t
  obtain ⟨g, hg⟩ := transitiveFromPoint t
  calc
    f t = f (actT g t0) := by rw [hg]
    _ = actU g (f t0) := fEquivariant g t0
    _ = actU g (h t0) := by rw [agreeAtPoint]
    _ = h (actT g t0) := by symm; exact hEquivariant g t0
    _ = h t := by rw [hg]

end FINPrograms295To308
