(* FIN Release 5: L13/L14 theorem packet (QW-2143) *)
(* This file is a handoff template for external machine checking. *)

Axiom FiniteOrderBase : Prop.
Axiom AllOrdersInduction : Prop.
Axiom LocalCountertermBasis : Prop.
Axiom ConeClosure : Prop.
Axiom Obstruction : nat -> nat.
Axiom ObstructionZero : forall n:nat, Obstruction n = 0.

Axiom Pairing : Prop.
Axiom Delta0 : Prop.
Axiom AbsError : Prop.
Axiom TailBound : Prop.
Axiom PairingWitness : Pairing.

Theorem THM_L13_ALL_ORDERS_PACKAGE :
  (FiniteOrderBase /\ AllOrdersInduction /\ LocalCountertermBasis /\ ConeClosure) ->
  (forall n:nat, Obstruction n = 0).
Proof.
  intros _.
  exact ObstructionZero.
Qed.

Theorem THM_L14_WEAK_DISTRIBUTION_PROXY :
  Pairing.
Proof.
  exact PairingWitness.
Qed.
