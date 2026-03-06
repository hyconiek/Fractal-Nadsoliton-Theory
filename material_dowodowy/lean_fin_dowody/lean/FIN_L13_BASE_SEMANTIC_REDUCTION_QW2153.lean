-- FIN Release 5: L13 base semantic reduction (QW-2153)
axiom Obstruction : Nat -> Nat

-- P2 obligation interpreted at theorem level:
axiom P2_finite_order_no_obstruction_n_le_4 :
  ∀ n : Nat, n <= 4 -> Obstruction n = 0

-- Remaining foundational bridge:
axiom P4_inductive_extension_rule : Prop
axiom step_from_P4 :
  P4_inductive_extension_rule -> (∀ n : Nat, Obstruction n = 0 -> Obstruction (Nat.succ n) = 0)

theorem base_from_P2_derived : Obstruction 0 = 0 := by
  have hle : (0 : Nat) <= 4 := by decide
  exact P2_finite_order_no_obstruction_n_le_4 0 hle

theorem THM_L13_STEP_ONLY_BRIDGE :
  P4_inductive_extension_rule ->
  (∀ n : Nat, Obstruction n = 0) := by
  intro h4
  have h0 : Obstruction 0 = 0 := base_from_P2_derived
  have hs : ∀ n : Nat, Obstruction n = 0 -> Obstruction (Nat.succ n) = 0 := step_from_P4 h4
  intro n
  induction n with
  | zero =>
    exact h0
  | succ k ih =>
    exact hs k ih
