-- FIN Release 5: L14 action-origin witness bundle (QW-2160)

theorem l14_action_origin_witness_bundle
  {c1 c2 c3 : Prop}
  (h1 : c1) (h2 : c2) (h3 : c3) :
  c1 ∧ c2 ∧ c3 := by
  exact And.intro h1 (And.intro h2 h3)
