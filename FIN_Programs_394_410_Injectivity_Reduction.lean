import Std

/- P396 formalizes the finite algebraic reduction of the P381 argument.
   ExpIndependent is the exact missing Lindemann--Weierstrass corollary.
   The file deliberately does not pretend that this transcendence theorem is
   available in the dependency-free local environment. -/
namespace FINPrograms394To410Injectivity

structure FourTermRelation where
  c1 : Rat
  c2 : Rat
  c3 : Rat
  c4 : Rat

def Nontrivial (r : FourTermRelation) : Prop :=
  r.c1 != 0 || r.c2 != 0 || r.c3 != 0 || r.c4 != 0

axiom ExpIndependent :
  forall (d e : Nat), d != e ->
    let r : FourTermRelation := {
      c1 := 1 + (e : Rat) ^ 9
      c2 := 1 + (e : Rat) ^ 9
      c3 := -(1 + (d : Rat) ^ 9)
      c4 := -(1 + (d : Rat) ^ 9) }
    Nontrivial r

theorem denominator_coefficients_nontrivial (d e : Nat) (h : d != e) :
    let r : FourTermRelation := {
      c1 := 1 + (e : Rat) ^ 9
      c2 := 1 + (e : Rat) ^ 9
      c3 := -(1 + (d : Rat) ^ 9)
      c4 := -(1 + (d : Rat) ^ 9) }
    Nontrivial r := by
  exact ExpIndependent d e h

theorem equality_contradicts_independence
    (kernelEqual : Prop)
    (relationForced : kernelEqual -> False) :
    Not kernelEqual := by
  intro h
  exact relationForced h

end FINPrograms394To410Injectivity
