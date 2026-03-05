-- FIN Release 5: L13 induction bridge decomposition (QW-2151)
axiom Obstruction : Nat -> Nat

-- Obligation-level carriers from prior strict chain:
axiom P2_finite_order_no_obstruction_n_le_4 : Prop
axiom P4_inductive_extension_rule : Prop

-- Remaining foundational assumptions decomposed into base and step:
axiom base_from_P2 : P2_finite_order_no_obstruction_n_le_4 -> Obstruction 0 = 0
axiom step_from_P4 :
  P4_inductive_extension_rule -> (∀ n : Nat, Obstruction n = 0 -> Obstruction (Nat.succ n) = 0)

theorem all_zero_from_base_step :
  (Obstruction 0 = 0) ->
  (∀ n : Nat, Obstruction n = 0 -> Obstruction (Nat.succ n) = 0) ->
  (∀ n : Nat, Obstruction n = 0) := by
  intro h0 hs
  intro n
  induction n with
  | zero =>
    exact h0
  | succ k ih =>
    exact hs k ih

theorem THM_L13_INDUCTION_REDUCED :
  P2_finite_order_no_obstruction_n_le_4 ->
  P4_inductive_extension_rule ->
  (∀ n : Nat, Obstruction n = 0) := by
  intro h2 h4
  have h0 : Obstruction 0 = 0 := base_from_P2 h2
  have hs : ∀ n : Nat, Obstruction n = 0 -> Obstruction (Nat.succ n) = 0 := step_from_P4 h4
  exact all_zero_from_base_step h0 hs
