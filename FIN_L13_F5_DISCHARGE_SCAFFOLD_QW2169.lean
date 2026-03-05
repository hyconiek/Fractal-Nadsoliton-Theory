-- FIN Release 5: L13 F5 discharge scaffold (QW-2169)

theorem l13_f5_discharge_scaffold
  {F5a F5b : Prop}
  (ha : F5a) (hb : F5b) :
  F5a ∧ F5b := by
  exact And.intro ha hb

theorem l13_f5_composition
  {F5a F5b F5 : Prop}
  (hcomp : F5a -> F5b -> F5)
  (ha : F5a) (hb : F5b) : F5 := by
  exact hcomp ha hb
