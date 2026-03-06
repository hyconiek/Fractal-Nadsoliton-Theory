-- FIN Release 5: L13 F5b terminal bound reduction (QW-2171)

theorem l13_f5b_conditional_bundle
  {A1 A2 A3 : Prop}
  (h1 : A1) (h2 : A2) (h3 : A3) :
  A1 ∧ A2 ∧ A3 := by
  exact And.intro h1 (And.intro h2 h3)

theorem l13_f5b_conditional_implies_f5b
  {A1 A2 A3 F5b : Prop}
  (hcond : A1 -> A2 -> A3 -> F5b)
  (h1 : A1) (h2 : A2) (h3 : A3) : F5b := by
  exact hcond h1 h2 h3
