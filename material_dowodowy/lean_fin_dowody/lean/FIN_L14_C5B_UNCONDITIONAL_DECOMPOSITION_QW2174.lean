-- FIN Release 5: L14 C5b unconditional decomposition (QW-2174)

theorem l14_c5b_v1_v2_bundle
  {V1 V2 : Prop}
  (h1 : V1) (h2 : V2) : V1 ∧ V2 := by
  exact And.intro h1 h2

theorem l14_c5b_from_v1_v2
  {V1 V2 C5b : Prop}
  (hcomp : V1 -> V2 -> C5b)
  (h1 : V1) (h2 : V2) : C5b := by
  exact hcomp h1 h2
