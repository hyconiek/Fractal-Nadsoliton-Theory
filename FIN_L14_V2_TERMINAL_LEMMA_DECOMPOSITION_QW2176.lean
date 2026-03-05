-- FIN Release 5: L14 V2 terminal lemma decomposition (QW-2176)

theorem l14_v2a_v2b_bundle
  {V2a V2b : Prop}
  (ha : V2a) (hb : V2b) : V2a ∧ V2b := by
  exact And.intro ha hb

theorem l14_v2_from_v2a_v2b
  {V2a V2b V2 : Prop}
  (hcomp : V2a -> V2b -> V2)
  (ha : V2a) (hb : V2b) : V2 := by
  exact hcomp ha hb
