-- FIN Release 5: L13 exhaustive canonical EoM gate (QW-2165)

theorem l13_exhaustive_canonical_eom_bundle
  {a b c d e f g : Prop}
  (ha : a) (hb : b) (hc : c) (hd : d) (he : e) (hf : f) (hg : g) :
  a ∧ b ∧ c ∧ d ∧ e ∧ f ∧ g := by
  exact And.intro ha (And.intro hb (And.intro hc (And.intro hd (And.intro he (And.intro hf hg)))))
