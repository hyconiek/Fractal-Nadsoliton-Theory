-- FIN Release 5: L13 U2 terminal lemma decomposition (QW-2175)

theorem l13_u2a_u2b_bundle
  {U2a U2b : Prop}
  (ha : U2a) (hb : U2b) : U2a ∧ U2b := by
  exact And.intro ha hb

theorem l13_u2_from_u2a_u2b
  {U2a U2b U2 : Prop}
  (hcomp : U2a -> U2b -> U2)
  (ha : U2a) (hb : U2b) : U2 := by
  exact hcomp ha hb
