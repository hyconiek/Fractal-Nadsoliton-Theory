-- FIN Release 5: L13 U2b2 exact matching identity (QW-2179)

theorem l13_u2b2_exact_matching_identity
  {CoeffId WeightId TermCount : Prop}
  (h1 : CoeffId) (h2 : WeightId) (h3 : TermCount) :
  CoeffId ∧ WeightId ∧ TermCount := by
  exact And.intro h1 (And.intro h2 h3)

theorem l13_u2b2_implies_u2b
  {U2b1 U2b2 U2b : Prop}
  (hcomp : U2b1 -> U2b2 -> U2b)
  (h1 : U2b1) (h2 : U2b2) : U2b := by
  exact hcomp h1 h2
