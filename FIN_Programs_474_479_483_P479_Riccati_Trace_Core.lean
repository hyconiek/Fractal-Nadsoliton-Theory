import Std

/-!
P479: dependency-free algebraic core of the O181 argument.

The carrier and operations are deliberately abstract.  Positivity, square
roots, spectral absolute values, and finite-matrix trace theory remain named
premises; this file proves only the exact implication structure actually used
by the FIN checkpoint.  It exports no selector, dimensions, physical units,
legacy-role transfer, L_total, Standard Model, GR, or laboratory claim.
-/

namespace FINPrograms474To483

theorem riccati_support_bridge
    {Matrix : Type}
    (riccati positiveSupport : Matrix -> Matrix -> Matrix -> Prop)
    (N X Delta : Matrix)
    (bridge : forall A B C, riccati A B C -> positiveSupport A B C)
    (hRiccati : riccati N X Delta) :
    positiveSupport N X Delta :=
  bridge N X Delta hRiccati

theorem trace_telescope_attainment
    {Primal Dual : Type}
    (primalValue : Primal -> Rat)
    (dualValue : Dual -> Rat)
    (feasiblePrimal : Primal -> Prop)
    (feasibleDual : Dual -> Prop)
    (contact : Primal -> Dual -> Prop)
    (weakDuality : forall p d,
      feasiblePrimal p -> feasibleDual d -> primalValue p <= dualValue d)
    (contactEquality : forall p d, contact p d -> primalValue p = dualValue d)
    (p : Primal) (d : Dual)
    (hp : feasiblePrimal p) (hd : feasibleDual d) (hc : contact p d) :
    primalValue p = dualValue d := by
  exact contactEquality p d hc

theorem exact_attainment_blocks_strict_improvement
    {Primal Dual : Type}
    (primalValue : Primal -> Rat)
    (dualValue : Dual -> Rat)
    (feasiblePrimal : Primal -> Prop)
    (feasibleDual : Dual -> Prop)
    (weakDuality : forall p d,
      feasiblePrimal p -> feasibleDual d -> primalValue p <= dualValue d)
    (witness : Primal) (certificate : Dual)
    (hw : feasiblePrimal witness) (hc : feasibleDual certificate)
    (attained : primalValue witness = dualValue certificate) :
    forall candidate, feasiblePrimal candidate ->
      primalValue candidate <= primalValue witness := by
  intro candidate hcandidate
  rw [attained]
  exact weakDuality candidate certificate hcandidate hc

theorem local_root_uniqueness_does_not_assert_global_optimizer_uniqueness
    (LocalRootUnique GlobalOptimizerUnique : Prop)
    (globalProof : GlobalOptimizerUnique) :
    GlobalOptimizerUnique :=
  globalProof

end FINPrograms474To483
