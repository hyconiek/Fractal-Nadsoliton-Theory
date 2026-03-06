-- FIN Release 5.1: QW-2365 L5 kernel-identity-integration discharge attempt
-- Scope: derive identity-integration theorem from identity-unification theorem.

variable (FINActionComplete ConstructiveNonPerturbativeScheme PositivityToReconstruction : Prop)

theorem QFT_KernelIdentityIntegrationToPositivity_Theorem :
  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by
  exact QFT_KernelIdentityUnificationToPositivity_Theorem

