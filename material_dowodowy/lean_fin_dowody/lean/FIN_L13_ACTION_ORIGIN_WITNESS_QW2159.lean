-- FIN Release 5: L13 action-origin witness bundle (QW-2159)

theorem l13_action_origin_witness_bundle
  {s1 s2 s3 s4 : Prop}
  (h1 : s1) (h2 : s2) (h3 : s3) (h4 : s4) :
  s1 ∧ s2 ∧ s3 ∧ s4 := by
  exact And.intro h1 (And.intro h2 (And.intro h3 h4))
