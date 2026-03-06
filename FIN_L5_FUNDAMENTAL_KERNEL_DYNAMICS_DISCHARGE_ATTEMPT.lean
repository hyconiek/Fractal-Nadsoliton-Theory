-- FIN Release 5.1: QW-2302 L5 fundamental-kernel-dynamics discharge attempt
-- Scope: derive kernel dynamics theorem from kernel-operator-closure theorem.

variable (FINActionComplete ConstructiveNonPerturbativeScheme PositivityToReconstruction : Prop)

theorem QFT_FundamentalKernelDynamicsToPositivity_Theorem :
  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by
  exact QFT_KernelOperatorClosureToPositivity_Theorem

