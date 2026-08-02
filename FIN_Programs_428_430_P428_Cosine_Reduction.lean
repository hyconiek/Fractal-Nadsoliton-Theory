import Std

/-!
P428: dependency-free reduction of the twelve FIN cosine enclosures.

The P411 interface used `Rat -> Rat`, which cannot literally denote the
standard real cosine at nonzero rational arguments.  This file repairs that
typing error by representing a real-like value only through rational lower
and upper predicates.  The local Lean/Std stack proves every FIN-specific
rational side condition.  A standard real-analysis theorem identifying the
cosine power series with the usual real cosine remains an external provider.
-/

namespace FINPrograms428To430Cosine

def factorial : Nat -> Nat
  | 0 => 1
  | n + 1 => (n + 1) * factorial n

def signedTerm (x : Rat) (k : Nat) : Rat :=
  let magnitude := x ^ (2*k) / (factorial (2*k) : Rat)
  if k % 2 = 0 then magnitude else -magnitude

def taylorSum (x : Rat) (n : Nat) : Rat :=
  (List.range (n+1)).foldl (fun total k => total + signedTerm x k) 0

def lower (x : Rat) : Rat := taylorSum x 21
def upper (x : Rat) : Rat := taylorSum x 20

def angles : List Rat := [
  ((650 : Rat) / 4000),
  ((1393 : Rat) / 4000),
  ((2136 : Rat) / 4000),
  ((2879 : Rat) / 4000),
  ((3622 : Rat) / 4000),
  ((4365 : Rat) / 4000),
  ((5108 : Rat) / 4000),
  ((5851 : Rat) / 4000),
  ((6594 : Rat) / 4000),
  ((7337 : Rat) / 4000),
  ((8080 : Rat) / 4000),
  ((8823 : Rat) / 4000)
]

def rationalSideConditions : Bool :=
  angles.all (fun x => decide (
    0 <= x && x*x < 12 &&
    lower x < upper x &&
    upper x - lower x < (1 : Rat) / 10^30))

theorem all_twelve_rational_side_conditions :
    rationalSideConditions = true := by
  native_decide

/- A rational-cut interface: unlike `Rat -> Rat`, this can type rational
   enclosures of values that are not rational numbers. -/
structure CutValue where
  lowerBound : Rat -> Prop
  upperBound : Rat -> Prop
  compatible : forall l u, lowerBound l -> upperBound u -> l <= u

abbrev CosineCutProvider := Rat -> CutValue

def StandardAlternatingCosineProvider
    (cosine : CosineCutProvider) : Prop :=
  forall x, 0 <= x -> x*x < 12 ->
    (cosine x).lowerBound (lower x) ∧
    (cosine x).upperBound (upper x)

def TwelveFINCosineBounds
    (cosine : CosineCutProvider) : Prop :=
  forall x, x ∈ angles ->
    (cosine x).lowerBound (lower x) ∧
    (cosine x).upperBound (upper x)

/- This theorem exposes the exact remaining dependency.  It is intentionally
   not an axiom claiming that a locally unavailable real cosine library exists. -/
theorem twelve_bounds_from_standard_provider
    (cosine : CosineCutProvider)
    (provider : StandardAlternatingCosineProvider cosine) :
    TwelveFINCosineBounds cosine := by
  intro x hx
  simp [angles] at hx
  rcases hx with rfl | rfl | rfl | rfl | rfl | rfl |
    rfl | rfl | rfl | rfl | rfl | rfl
  all_goals exact provider _ (by native_decide) (by native_decide)

end FINPrograms428To430Cosine
