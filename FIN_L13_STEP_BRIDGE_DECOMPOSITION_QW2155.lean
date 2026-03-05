-- FIN Release 5: L13 step bridge decomposition (QW-2155)
axiom Obstruction : Nat -> Nat
axiom P4_inductive_extension_rule : Prop

-- Step bridge decomposed into explicit sub-obligations:
axiom step_s1_local_counterterm_lift : Prop
axiom step_s2_weighted_remainder_contractive : Prop
axiom step_s3_distribution_split_stable : Prop
axiom step_s4_obstruction_projection_zero : Prop

def StepBridgeBundle : Prop :=
  step_s1_local_counterterm_lift ∧
  step_s2_weighted_remainder_contractive ∧
  step_s3_distribution_split_stable ∧
  step_s4_obstruction_projection_zero

axiom p4_implies_step_bundle : P4_inductive_extension_rule -> StepBridgeBundle
axiom step_bundle_implies_step_rule :
  StepBridgeBundle -> (∀ n : Nat, Obstruction n = 0 -> Obstruction (Nat.succ n) = 0)

theorem step_rule_from_p4_decomposed :
  P4_inductive_extension_rule -> (∀ n : Nat, Obstruction n = 0 -> Obstruction (Nat.succ n) = 0) := by
  intro h4
  have hb : StepBridgeBundle := p4_implies_step_bundle h4
  exact step_bundle_implies_step_rule hb
