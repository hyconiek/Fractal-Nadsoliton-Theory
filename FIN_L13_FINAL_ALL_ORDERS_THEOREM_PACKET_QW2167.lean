-- FIN Release 5: L13 final all-orders theorem packet (QW-2167)

theorem l13_final_all_orders_packet
  {F1 F2 F3 F4 F5 : Prop}
  (h1 : F1) (h2 : F2) (h3 : F3) (h4 : F4)
  (h5 : F5) :
  F1 ∧ F2 ∧ F3 ∧ F4 ∧ F5 := by
  exact And.intro h1 (And.intro h2 (And.intro h3 (And.intro h4 h5)))
