import Std

/-!
FIN Programs P281--P294: dependency-free exact Schur certificate.

This file is deliberately limited to Lean core and exact rational arithmetic.
It machine-checks:

1. identity reduction on a retained scalar block;
2. equality of nested and direct Schur elimination for a frozen rational
   positive 3 x 3 witness;
3. positivity of every pivot and of the final reduced scalar in that witness;
4. equality of the two reduction paths as an exact rational statement.

The file is not a Mathlib formalization of the general positive-matrix Schur
theorem.  That larger theorem remains a separate recommended program.
-/

structure Symmetric3 where
  a : Rat
  b : Rat
  c : Rat
  d : Rat
  e : Rat
  f : Rat
deriving Repr, DecidableEq

def schurIdentity (x : Rat) : Rat := x

def eliminateThird (m : Symmetric3) : Rat × Rat × Rat :=
  (m.a - m.c * m.c / m.f,
   m.b - m.c * m.e / m.f,
   m.d - m.e * m.e / m.f)

def nestedSchurToFirst (m : Symmetric3) : Rat :=
  let reduced := eliminateThird m
  reduced.1 - reduced.2.1 * reduced.2.1 / reduced.2.2

def directSchurToFirst (m : Symmetric3) : Rat :=
  m.a -
    (m.f * m.b * m.b - 2 * m.e * m.b * m.c + m.d * m.c * m.c) /
      (m.d * m.f - m.e * m.e)

def witness : Symmetric3 where
  a := 3
  b := -1
  c := 0
  d := 3
  e := -1
  f := 2

theorem identity_reduction_exact :
    schurIdentity witness.a = witness.a := by
  rfl

theorem third_pivot_positive :
    0 < witness.f := by
  native_decide

theorem hidden_block_determinant_positive :
    0 < witness.d * witness.f - witness.e * witness.e := by
  native_decide

theorem intermediate_pivot_positive :
    0 < (eliminateThird witness).2.2 := by
  native_decide

theorem nested_value_exact :
    nestedSchurToFirst witness = 13 / 5 := by
  native_decide

theorem direct_value_exact :
    directSchurToFirst witness = 13 / 5 := by
  native_decide

theorem nested_direct_composition_exact :
    nestedSchurToFirst witness = directSchurToFirst witness := by
  native_decide

theorem reduced_value_positive :
    0 < nestedSchurToFirst witness := by
  native_decide

#eval nestedSchurToFirst witness
#eval directSchurToFirst witness
