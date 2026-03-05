-- FIN Release 5: L14 continuum sub-obligation grounding (QW-2158)

theorem continuum_bundle_from_grounded_conditions
  {b11 b12 b13 b21 b22 b23 b31 b32 b33 : Prop}
  (h11 : b11) (h12 : b12) (h13 : b13)
  (h21 : b21) (h22 : b22) (h23 : b23)
  (h31 : b31) (h32 : b32) (h33 : b33) :
  (b11 ∧ b12 ∧ b13) ∧
  (b21 ∧ b22 ∧ b23) ∧
  (b31 ∧ b32 ∧ b33) := by
  refine And.intro (And.intro h11 (And.intro h12 h13)) ?_
  refine And.intro (And.intro h21 (And.intro h22 h23)) ?_
  exact And.intro h31 (And.intro h32 h33)
