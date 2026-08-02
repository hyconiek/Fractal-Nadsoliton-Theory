import Std

/- Exact rational reflection of the Taylor polynomials used to enclose
   cos((743*d+650)/4000), d=0,...,11.  The analytic alternating-remainder
   theorem is a separately named external dependency. -/
namespace FINPrograms411To427Taylor

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

def rationalCertificate : Bool :=
  angles.all (fun x => decide (lower x < upper x && upper x - lower x < (1 : Rat) / 10^30))

theorem all_twelve_rational_taylor_checks : rationalCertificate = true := by
  native_decide

def AnalyticCosineBridge (cosine : Rat -> Rat) : Prop :=
  forall x, x ∈ angles -> lower x <= cosine x && cosine x <= upper x

theorem interval_use_requires_bridge
    (cosine : Rat -> Rat) (bridge : AnalyticCosineBridge cosine) :
    forall x, x ∈ angles -> lower x <= cosine x && cosine x <= upper x :=
  bridge

end FINPrograms411To427Taylor
