-- FIN Release 5: L14 C5 discharge scaffold (QW-2170)

theorem l14_c5_discharge_scaffold
  {C5a C5b : Prop}
  (ha : C5a) (hb : C5b) :
  C5a ∧ C5b := by
  exact And.intro ha hb

theorem l14_c5_composition
  {C5a C5b C5 : Prop}
  (hcomp : C5a -> C5b -> C5)
  (ha : C5a) (hb : C5b) : C5 := by
  exact hcomp ha hb
