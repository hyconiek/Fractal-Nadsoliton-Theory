import Std

/-!
FIN Programs P323--P336: dependency-free formal core.

The statements below formalize only the structural implications used by
Release 10.29.  They do not derive numerical FIN parameters, a physical
sector, a dimensional unit, an apparatus, or a legacy-to-strict bridge.
-/

namespace FINPrograms323To336

theorem source_independent_completion_nonidentification
    {Legacy Target : Type}
    (complete : Legacy → Target)
    (legacy : Legacy)
    (target₁ target₂ : Target)
    (different : target₁ ≠ target₂) :
    ¬ (complete legacy = target₁ ∧ complete legacy = target₂) := by
  intro both
  apply different
  exact both.1.symm.trans both.2

theorem dual_lower_bound
    {Measure Polynomial Scalar : Type}
    (cost : Measure → Scalar)
    (value : Polynomial → Scalar)
    (le : Scalar → Scalar → Prop)
    (feasibleMeasure : Measure → Prop)
    (feasiblePolynomial : Polynomial → Prop)
    (weakDuality :
      ∀ μ p, feasibleMeasure μ → feasiblePolynomial p →
        le (value p) (cost μ))
    (μ : Measure)
    (p : Polynomial)
    (hμ : feasibleMeasure μ)
    (hp : feasiblePolynomial p) :
    le (value p) (cost μ) :=
  weakDuality μ p hμ hp

theorem matching_primal_dual_is_optimal
    {Measure Polynomial Scalar : Type}
    (cost : Measure → Scalar)
    (value : Polynomial → Scalar)
    (le : Scalar → Scalar → Prop)
    (feasibleMeasure : Measure → Prop)
    (feasiblePolynomial : Polynomial → Prop)
    (weakDuality :
      ∀ μ p, feasibleMeasure μ → feasiblePolynomial p →
        le (value p) (cost μ))
    (μstar : Measure)
    (pstar : Polynomial)
    (hμ : feasibleMeasure μstar)
    (hp : feasiblePolynomial pstar)
    (matchingValue : value pstar = cost μstar) :
    (∀ μ, feasibleMeasure μ → le (cost μstar) (cost μ)) ∧
    (∀ p, feasiblePolynomial p → le (value p) (value pstar)) := by
  constructor
  · intro μ h
    rw [← matchingValue]
    exact weakDuality μ pstar h hp
  · intro p h
    rw [matchingValue]
    exact weakDuality μstar p hμ h

theorem two_sector_swap_has_no_fixed_point :
    ¬ ∃ sector : Bool, (!sector) = sector := by
  intro fixed
  obtain ⟨sector, equality⟩ := fixed
  cases sector <;> cases equality

abbrev ResourceVector := Fin 5 → Bool

def allExcept (missing : Fin 5) : ResourceVector :=
  fun index => decide (index ≠ missing)

theorem independence_cube_witness (missing : Fin 5) :
    (allExcept missing missing = false) ∧
    (∀ index, index ≠ missing → allExcept missing index = true) := by
  constructor
  · simp [allExcept]
  · intro index different
    simp [allExcept, different]

theorem no_single_resource_follows_from_the_other_four
    (missing : Fin 5) :
    ∃ resources : ResourceVector,
      resources missing = false ∧
      ∀ index, index ≠ missing → resources index = true := by
  exact ⟨allExcept missing, independence_cube_witness missing⟩

end FINPrograms323To336
