-- FIN Release 5: L14 final continuum theorem packet (QW-2168)

theorem l14_final_continuum_packet
  {C1 C2 C3 C4 C5 : Prop}
  (h1 : C1) (h2 : C2) (h3 : C3) (h4 : C4)
  (h5 : C5) :
  C1 ∧ C2 ∧ C3 ∧ C4 ∧ C5 := by
  exact And.intro h1 (And.intro h2 (And.intro h3 (And.intro h4 h5)))
