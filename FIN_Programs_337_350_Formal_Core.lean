import Std

/-!
FIN Programs P337--P350: dependency-free structural core.

These theorems encode only the abstract implications used in Release 10.30.
They do not derive FIN numerical parameters, an internal selector, physical
units, apparatus data, electroweak dynamics, a legacy-role transfer theorem,
L_total, SM/GR, or a Theory of Everything.
-/

namespace FINPrograms337To350

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

theorem source_independent_target_nonidentification
    {Input Target : Type}
    (complete : Input → Target)
    (input : Input)
    (target₁ target₂ : Target)
    (different : target₁ ≠ target₂) :
    ¬ (complete input = target₁ ∧ complete input = target₂) := by
  intro both
  apply different
  exact both.1.symm.trans both.2

theorem commuting_function_preserves_commutator_zero
    {Operator : Type}
    (compose : Operator → Operator → Operator)
    (source mapped : Operator)
    (commutes : compose source mapped = compose mapped source) :
    compose source mapped = compose mapped source :=
  commutes

theorem finite_observation_nonidentification
    {Clock Observation : Type}
    (observe : Clock → Observation)
    (clock₁ clock₂ : Clock)
    (sameRecord : observe clock₁ = observe clock₂)
    (different : clock₁ ≠ clock₂) :
    ∃ first second,
      first ≠ second ∧ observe first = observe second :=
  ⟨clock₁, clock₂, different, sameRecord⟩

inductive BridgeResource where
  | signedPath
  | nontorsionPhase
  | crossLaw
  | pointing
  | dimensionalScale
  deriving DecidableEq

abbrev ResourceVector := BridgeResource → Bool

def allExcept (missing : BridgeResource) : ResourceVector :=
  fun resource => decide (resource ≠ missing)

theorem typed_single_omission_witness (missing : BridgeResource) :
    (allExcept missing missing = false) ∧
    (∀ resource, resource ≠ missing → allExcept missing resource = true) := by
  constructor
  · simp [allExcept]
  · intro resource different
    simp [allExcept, different]

theorem conditional_result_keeps_assumption
    {Assumption Conclusion : Prop}
    (derivation : Assumption → Conclusion)
    (assumption : Assumption) :
    Conclusion :=
  derivation assumption

end FINPrograms337To350
