-- FIN Release 5: L14 V2b action-bridge spec (QW-2178)

theorem l14_v2b1_v2b2_bundle
  {V2b1 V2b2 : Prop}
  (h1 : V2b1) (h2 : V2b2) : V2b1 ∧ V2b2 := by
  exact And.intro h1 h2

theorem l14_v2b_from_v2b1_v2b2
  {V2b1 V2b2 V2b : Prop}
  (hcomp : V2b1 -> V2b2 -> V2b)
  (h1 : V2b1) (h2 : V2b2) : V2b := by
  exact hcomp h1 h2
