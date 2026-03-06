-- FIN Release 5.1: QW-2362 L5 kernel-identity-consolidation discharge attempt
-- Scope: derive identity-consolidation theorem from identity-integration theorem.

variable (FINActionComplete ConstructiveNonPerturbativeScheme PositivityToReconstruction : Prop)

theorem QFT_KernelIdentityConsolidationToPositivity_Theorem :
  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by
  exact QFT_KernelIdentityIntegrationToPositivity_Theorem

