import Std

/-!
FIN Programs P309--P322: dependency-free formal core.

The file checks three structural statements used in this release:

* P310: a completion-of-the-square identity plus a nonnegative remainder
  proves the Schur minimizer;
* P314/P322: a finite/torsion source cannot hit a declared nontorsion target
  under a torsion-preserving map, and the two-sector swap has no invariant
  point;
* P322: every observable that factors through a kernel is constant on the
  fibers of that kernel.

The hypotheses do not derive a square decomposition, a phase resource, a
sector point, a legacy-to-strict bridge, or a physical role.
-/

namespace FINPrograms309To322

theorem completion_square_minimizer
    {X Y V : Type}
    (add : V → V → V)
    (zero : V)
    (addZero : ∀ a : V, add a zero = a)
    (le : V → V → Prop)
    (addMonoLeft : ∀ a b c : V, le b c → le (add a b) (add a c))
    (F : X → Y → V)
    (reduced : X → V)
    (remainder : X → Y → V)
    (yStar : X → Y)
    (decomposition :
      ∀ x y, F x y = add (reduced x) (remainder x y))
    (remainderNonnegative :
      ∀ x y, le zero (remainder x y))
    (remainderVanishes :
      ∀ x, remainder x (yStar x) = zero) :
    ∀ x y, le (F x (yStar x)) (F x y) := by
  intro x y
  rw [decomposition x (yStar x), remainderVanishes x, addZero]
  rw [decomposition x y]
  have h := remainderNonnegative x y
  simpa only [addZero] using
    addMonoLeft (reduced x) zero (remainder x y) h

theorem torsion_preserving_map_misses_nontorsion_target
    {S T : Type}
    (torsionS : S → Prop)
    (torsionT : T → Prop)
    (f : S → T)
    (sourceTorsion : ∀ s, torsionS s)
    (preservesTorsion : ∀ s, torsionS s → torsionT (f s))
    (target : T)
    (targetNontorsion : ¬ torsionT target) :
    ∀ s, f s ≠ target := by
  intro s equality
  apply targetNontorsion
  rw [← equality]
  exact preservesTorsion s (sourceTorsion s)

theorem bool_swap_has_no_invariant_point :
    ¬ ∃ b : Bool, (!b) = b := by
  intro fixed
  obtain ⟨b, hb⟩ := fixed
  cases b <;> cases hb

theorem factorized_role_constant_on_kernel_fibers
    {P K O : Type}
    (kernel : P → K)
    (roleOnKernel : K → O)
    (role : P → O)
    (factorization : ∀ p, role p = roleOnKernel (kernel p)) :
    ∀ p q, kernel p = kernel q → role p = role q := by
  intro p q kernelEquality
  rw [factorization p, factorization q, kernelEquality]

end FINPrograms309To322
