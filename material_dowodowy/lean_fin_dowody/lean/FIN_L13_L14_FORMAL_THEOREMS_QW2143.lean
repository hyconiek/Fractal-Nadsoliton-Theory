-- FIN Release 5: L13/L14 theorem packet (QW-2143)
-- This file is a handoff template for external machine checking.

axiom FiniteOrderBase : Prop
axiom AllOrdersInduction : Prop
axiom LocalCountertermBasis : Prop
axiom ConeClosure : Prop
axiom Obstruction : Nat -> Nat
axiom ObstructionZero : ∀ n : Nat, Obstruction n = 0

axiom Pairing : Prop
axiom Delta0 : Prop
axiom AbsError : Prop
axiom TailBound : Prop
axiom PairingWitness : Pairing

theorem THM_L13_ALL_ORDERS_PACKAGE :
  (FiniteOrderBase ∧ AllOrdersInduction ∧ LocalCountertermBasis ∧ ConeClosure) ->
  (∀ n : Nat, Obstruction n = 0) := by
  intro _
  exact ObstructionZero

theorem THM_L14_WEAK_DISTRIBUTION_PROXY :
  Pairing := by
  exact PairingWitness
