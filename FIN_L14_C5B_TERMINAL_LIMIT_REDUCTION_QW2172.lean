-- FIN Release 5: L14 C5b terminal limit reduction (QW-2172)

theorem l14_c5b_conditional_bundle
  {B1 B2 B3 : Prop}
  (h1 : B1) (h2 : B2) (h3 : B3) :
  B1 ∧ B2 ∧ B3 := by
  exact And.intro h1 (And.intro h2 h3)

theorem l14_c5b_conditional_implies_c5b
  {B1 B2 B3 C5b : Prop}
  (hcond : B1 -> B2 -> B3 -> C5b)
  (h1 : B1) (h2 : B2) (h3 : B3) : C5b := by
  exact hcond h1 h2 h3
