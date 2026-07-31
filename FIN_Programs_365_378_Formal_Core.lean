import Std

/-!
Dependency-free structural core for FIN Programs P365--P378.

The exact numerical certificates are checked by separate arbitrary-precision
executables. This Lean file formalizes only the structural implications and
does not claim proof-assistant reflection of their large rational arithmetic.
-/

namespace FINPrograms365To378

abbrev Vec5 := Fin 5 → Nat

def addVec (left right : Vec5) : Vec5 :=
  fun index => left index + right index

def Free (source target : Vec5) : Prop :=
  ∀ index, target index ≤ source index

theorem free_composition
    (a b c : Vec5)
    (hab : Free a b)
    (hbc : Free b c) :
    Free a c := by
  intro index
  exact Nat.le_trans (hbc index) (hab index)

theorem coordinate_additivity
    (left right : Vec5) (index : Fin 5) :
    addVec left right index = left index + right index := rfl

theorem maximal_kernel_morphism_definition
    {Vertex Distance Value : Type}
    (distanceSource distanceTarget : Vertex → Vertex → Distance)
    (kernel : Distance → Value)
    (map : Vertex → Vertex)
    (preserves :
      ∀ x y,
        kernel (distanceTarget (map x) (map y)) =
        kernel (distanceSource x y)) :
    ∀ x y,
      kernel (distanceTarget (map x) (map y)) =
      kernel (distanceSource x y) :=
  preserves

theorem angle_upper_lower_collapse
    {Value : Type}
    (le : Value → Value → Prop)
    (adaptive parallel upper : Value)
    (parallelFeasible : le parallel adaptive)
    (adaptiveBound : le adaptive upper)
    (matching : parallel = upper) :
    le parallel adaptive ∧ le adaptive parallel := by
  constructor
  · exact parallelFeasible
  · rw [matching]
    exact adaptiveBound

theorem data_processing_instance
    {State : Type}
    (monotone : State → Nat)
    (channel : State → State)
    (contracts : ∀ state, monotone (channel state) ≤ monotone state) :
    ∀ state, monotone (channel state) ≤ monotone state :=
  contracts

theorem resource_grade_does_not_fix_realization
    {Grade Realization : Type}
    (grade : Realization → Grade)
    (left right : Realization)
    (sameGrade : grade left = grade right)
    (different : left ≠ right) :
    ∃ first second,
      grade first = grade second ∧ first ≠ second :=
  ⟨left, right, sameGrade, different⟩

end FINPrograms365To378
