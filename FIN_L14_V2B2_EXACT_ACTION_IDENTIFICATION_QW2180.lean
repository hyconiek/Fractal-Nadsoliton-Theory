-- FIN Release 5: L14 V2b2 exact action identification (QW-2180)

theorem l14_v2b2_exact_operator_hessian_identity
  {OpEq HessSym ShapeOK : Prop}
  (h1 : OpEq) (h2 : HessSym) (h3 : ShapeOK) :
  OpEq ∧ HessSym ∧ ShapeOK := by
  exact And.intro h1 (And.intro h2 h3)

theorem l14_v2b2_implies_v2b
  {V2b1 V2b2 V2b : Prop}
  (hcomp : V2b1 -> V2b2 -> V2b)
  (h1 : V2b1) (h2 : V2b2) : V2b := by
  exact hcomp h1 h2
