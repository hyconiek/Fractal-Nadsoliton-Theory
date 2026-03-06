-- FIN Release 5: L13 F5b unconditional decomposition (QW-2173)

theorem l13_f5b_u1_u2_bundle
  {U1 U2 : Prop}
  (h1 : U1) (h2 : U2) : U1 ∧ U2 := by
  exact And.intro h1 h2

theorem l13_f5b_from_u1_u2
  {U1 U2 F5b : Prop}
  (hcomp : U1 -> U2 -> F5b)
  (h1 : U1) (h2 : U2) : F5b := by
  exact hcomp h1 h2
