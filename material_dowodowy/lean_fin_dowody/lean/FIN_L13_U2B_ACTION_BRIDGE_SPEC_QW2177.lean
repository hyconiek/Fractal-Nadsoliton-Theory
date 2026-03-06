-- FIN Release 5: L13 U2b action-bridge spec (QW-2177)

theorem l13_u2b1_u2b2_bundle
  {U2b1 U2b2 : Prop}
  (h1 : U2b1) (h2 : U2b2) : U2b1 ∧ U2b2 := by
  exact And.intro h1 h2

theorem l13_u2b_from_u2b1_u2b2
  {U2b1 U2b2 U2b : Prop}
  (hcomp : U2b1 -> U2b2 -> U2b)
  (h1 : U2b1) (h2 : U2b2) : U2b := by
  exact hcomp h1 h2
