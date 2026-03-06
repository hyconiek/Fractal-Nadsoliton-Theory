-- FIN Release 5.1: QW-2374 L5 kernel-identity-totality discharge attempt
-- Scope: derive identity-totality theorem from identity-finality theorem.

variable (FINActionComplete ConstructiveNonPerturbativeScheme PositivityToReconstruction : Prop)

theorem QFT_KernelIdentityTotalityToPositivity_Theorem :
  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by
  exact QFT_KernelIdentityFinalityToPositivity_Theorem

