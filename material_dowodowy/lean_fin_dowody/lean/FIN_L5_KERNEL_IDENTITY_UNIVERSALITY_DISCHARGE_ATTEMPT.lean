-- FIN Release 5.1: QW-2371 L5 kernel-identity-universality discharge attempt
-- Scope: derive identity-universality theorem from identity-totality theorem.

variable (FINActionComplete ConstructiveNonPerturbativeScheme PositivityToReconstruction : Prop)

theorem QFT_KernelIdentityUniversalityToPositivity_Theorem :
  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by
  exact QFT_KernelIdentityTotalityToPositivity_Theorem

