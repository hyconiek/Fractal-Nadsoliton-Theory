-- FIN Release 5: dual L13/L14 terminal matching closure (QW-2181)

theorem l13_l14_dual_terminal_closure
  {L13Term L14Term : Prop}
  (h13 : L13Term) (h14 : L14Term) : L13Term ∧ L14Term := by
  exact And.intro h13 h14
