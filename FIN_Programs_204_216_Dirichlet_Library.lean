import Mathlib.Data.Matrix.Basic
import Mathlib.Data.Matrix.Notation
import Mathlib.LinearAlgebra.Matrix.PosDef
import Mathlib.Analysis.Complex.Exponential

/-!
FIN Programs 204--216: finite Dirichlet library.

This is the general Mathlib-backed source targeted by the P207 build
contract.  The release records the core Lean probe separately from the
compilation status of this full library.
-/

open scoped BigOperators Matrix

namespace FINSpectral

variable {n : Type} [Fintype n] [DecidableEq n]

def rowSum (W : Matrix n n ℝ) (x : n) : ℝ := ∑ y, W x y

def generator (s : ℝ) (W : Matrix n n ℝ) : Matrix n n ℝ :=
  s • (1 : Matrix n n ℝ) - W

theorem dirichlet_identity
    (W : Matrix n n ℝ)
    (hSymm : W.IsHermitian)
    (s : ℝ)
    (hRow : ∀ x, rowSum W x = s)
    (f : n → ℝ) :
    (∑ x, f x * (generator s W).mulVec f x)
      = (1 / 2 : ℝ) * ∑ x, ∑ y, W x y * (f x - f y)^2 := by
  classical
  simp only [generator, Matrix.sub_mulVec, Matrix.smul_mulVec, Matrix.one_mulVec]
  simp only [Matrix.mulVec, dotProduct]
  rw [Finset.mul_sum]
  simp_rw [Finset.mul_sum, hRow]
  have hswap :
      (∑ x, ∑ y, W x y * f y ^ 2) = ∑ x, s * f x ^ 2 := by
    calc
      (∑ x, ∑ y, W x y * f y ^ 2)
          = ∑ y, ∑ x, W x y * f y ^ 2 := by
              rw [Finset.sum_comm]
      _ = ∑ y, (∑ x, W x y) * f y ^ 2 := by
              apply Finset.sum_congr rfl
              intro y _
              rw [Finset.sum_mul]
      _ = ∑ y, s * f y ^ 2 := by
              apply Finset.sum_congr rfl
              intro y _
              have hcol : (∑ x, W x y) = s := by
                simpa [rowSum, hSymm.eq y] using hRow y
              rw [hcol]
  rw [show (1 / 2 : ℝ) * ∑ x, ∑ y, W x y * (f x - f y)^2
      = (1 / 2 : ℝ) *
        ((∑ x, s * f x^2) + (∑ x, s * f x^2)
          - 2 * ∑ x, f x * ∑ y, W x y * f y) by
        simp_rw [sub_sq, mul_add, mul_sub, Finset.sum_add_distrib,
          Finset.sum_sub_distrib]
        rw [hswap]
        ring]
  ring

theorem generator_nonnegative
    (W : Matrix n n ℝ)
    (hSymm : W.IsHermitian)
    (hNonneg : ∀ x y, 0 ≤ W x y)
    (s : ℝ)
    (hRow : ∀ x, rowSum W x = s)
    (f : n → ℝ) :
    0 ≤ ∑ x, f x * (generator s W).mulVec f x := by
  rw [dirichlet_identity W hSymm s hRow f]
  positivity

theorem constant_in_kernel
    (W : Matrix n n ℝ)
    (s c : ℝ)
    (hRow : ∀ x, rowSum W x = s) :
    (generator s W).mulVec (fun _ => c) = 0 := by
  funext x
  simp [generator, Matrix.mulVec, dotProduct, rowSum, hRow x,
    Finset.mul_sum]

end FINSpectral
