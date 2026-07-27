def energy2 (w f0 f1 : Nat) : Nat :=
  w * (f0 - f1) * (f0 - f1)

#eval energy2 3 5 2

example : energy2 3 5 2 = 27 := by decide
example : energy2 7 4 4 = 0 := by decide
