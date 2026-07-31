import Std

/-!
Structural formal core for FIN Programs P351--P364.

This dependency-free file proves only the abstract implications used in the
report.  The rational interval certificates are emitted by the executable
Python artifact and are not silently imported as Lean theorems.
-/

namespace FINPrograms351To364

def Vec5 := Fin 5 → Nat

def Free (r s : Vec5) : Prop := ∀ i, s i ≤ r i

theorem coordinate_monotones_complete (r s : Vec5) :
    Free r s ↔ ∀ i, s i ≤ r i := Iff.rfl

def addVec (r s : Vec5) : Vec5 := fun i => r i + s i

theorem additive_cancellation (r s c : Vec5) :
    Free (addVec r c) (addVec s c) ↔ Free r s := by
  constructor
  · intro h i
    exact Nat.le_of_add_le_add_right (h i)
  · intro h i
    exact Nat.add_le_add_right (h i) (c i)

theorem bounded_strategy_optimal
    {Value : Type}
    (le : Value → Value → Prop)
    (achieved upper : Value)
    (hAchieved : achieved = upper)
    (hBound : ∀ candidate : Value, le candidate upper) :
    ∀ candidate : Value, le candidate achieved := by
  intro candidate
  rw [hAchieved]
  exact hBound candidate

theorem radial_isomorphism_naturality
    {V : Type} (d : V → V → Nat) (k : Nat → ℝ)
    (p : V → V) (hIso : ∀ x y, d (p x) (p y) = d x y) :
    ∀ x y, k (d (p x) (p y)) = k (d x y) := by
  intro x y
  rw [hIso]

theorem minimax_design_from_gap_bound
    {Design : Type}
    (gap : Design → Nat)
    (equispaced : Design)
    (lowerBound : Nat)
    (allDesigns : ∀ design, lowerBound ≤ gap design)
    (attains : gap equispaced = lowerBound) :
    ∀ design, gap equispaced ≤ gap design := by
  intro design
  rw [attains]
  exact allDesigns design

structure TypedKernel where
  label : String

structure OperationalRecord where
  preparation : String
  channel : String
  instrument : String
  recordHash : String

def HasOperationalBridge (K : TypedKernel) : Prop :=
  Nonempty OperationalRecord

theorem mathematics_does_not_supply_record
    (K : TypedKernel) :
    HasOperationalBridge K → Nonempty OperationalRecord := by
  intro h
  exact h

end FINPrograms351To364
