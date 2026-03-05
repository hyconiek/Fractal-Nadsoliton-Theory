-- FIN Release 5: L13 full canonical action variational gate (QW-2163)

theorem l13_full_canonical_action_variational_bundle
  {a b c d e : Prop}
  (ha : a) (hb : b) (hc : c) (hd : d) (he : e) :
  a ∧ b ∧ c ∧ d ∧ e := by
  exact And.intro ha (And.intro hb (And.intro hc (And.intro hd he)))
