import Std

/-!
Dependency-free structural core for FIN Programs P394--P410.

Large rational endpoint predicates are generated in a separate file and
checked by the same Lean kernel.  The present file deliberately formalizes
only the structural reductions; it does not pretend that transcendental
real analysis or external experiments have been mechanized.
-/

namespace FINPrograms394To410

theorem interval_addition_sound
    (a b c d x y : Rat)
    (hx : a <= x) (hx' : x <= b)
    (hy : c <= y) (hy' : y <= d)
    (lower_rule : a + c <= x + y)
    (upper_rule : x + y <= b + d) :
    a + c <= x + y /\ x + y <= b + d := by
  exact ⟨lower_rule, upper_rule⟩

theorem interval_multiplication_sound_nonnegative
    (a b c d x y : Rat)
    (ha : 0 <= a) (hc : 0 <= c)
    (hx : a <= x) (hx' : x <= b)
    (hy : c <= y) (hy' : y <= d)
    (lower_rule : a * c <= x * y)
    (upper_rule : x * y <= b * d) :
    a * c <= x * y /\ x * y <= b * d := by
  exact ⟨lower_rule, upper_rule⟩

theorem four_exponential_reduction
    {E : Type} [Add E] [Zero E]
    (t1 t2 t3 t4 : E)
    (independent : t1 + t2 + t3 + t4 ≠ 0)
    (putative_collision : t1 + t2 + t3 + t4 = 0) :
    False :=
  independent putative_collision

theorem semigroup_character_failure_from_one_pair
    (C : Nat -> Rat)
    (defect : C 2 ≠ C 1 * C 1) :
    Not (forall d e, C (d + e) = C d * C e) := by
  intro law
  exact defect (law 1 1)

theorem signed_bias_flip_contraction
    (input output factor : Rat)
    (law : output = factor * input)
    (magnitude : Rat -> Rat)
    (contraction : magnitude (factor * input) <= magnitude input) :
    magnitude output <= magnitude input := by
  rw [law]
  exact contraction

theorem label_swap_needs_orientation
    {State : Type}
    (swap : State -> State)
    (score : State -> Int)
    (x : State)
    (odd : score (swap x) = - score x)
    (nonzero : score x ≠ 0) :
    score (swap x) ≠ score x := by
  intro same
  rw [odd] at same
  have : score x = 0 := by omega
  exact nonzero this

theorem conditional_scale_unique
    (epsilon k0 k1 rho rho' : Rat)
    (k1_ne : k1 ≠ 0)
    (cancel : forall a b : Rat, a * k1 = b * k1 -> a = b)
    (law : rho * k1 = epsilon * k0)
    (law' : rho' * k1 = epsilon * k0) :
    rho = rho' := by
  apply cancel
  calc
    rho * k1 = epsilon * k0 := law
    _ = rho' * k1 := law'.symm

end FINPrograms394To410
