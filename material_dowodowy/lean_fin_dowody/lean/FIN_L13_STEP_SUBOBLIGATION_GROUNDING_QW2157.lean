-- FIN Release 5: L13 step sub-obligation grounding (QW-2157)

theorem step_bundle_from_grounded_conditions
  {a11 a12 a13 a21 a22 a23 a31 a32 a33 a34 a41 a42 : Prop}
  (h11 : a11) (h12 : a12) (h13 : a13)
  (h21 : a21) (h22 : a22) (h23 : a23)
  (h31 : a31) (h32 : a32) (h33 : a33) (h34 : a34)
  (h41 : a41) (h42 : a42) :
  (a11 ∧ a12 ∧ a13) ∧
  (a21 ∧ a22 ∧ a23) ∧
  (a31 ∧ a32 ∧ a33 ∧ a34) ∧
  (a41 ∧ a42) := by
  refine And.intro (And.intro h11 (And.intro h12 h13)) ?_
  refine And.intro (And.intro h21 (And.intro h22 h23)) ?_
  refine And.intro (And.intro h31 (And.intro h32 (And.intro h33 h34))) ?_
  exact And.intro h41 h42
