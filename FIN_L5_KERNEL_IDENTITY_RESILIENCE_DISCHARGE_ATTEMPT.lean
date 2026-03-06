-- FIN Release 5.1: QW-2359 L5 kernel-identity-resilience discharge attempt
-- Scope: derive identity-resilience theorem from identity-consolidation theorem.

variable (FINActionComplete ConstructiveNonPerturbativeScheme PositivityToReconstruction : Prop)

theorem QFT_KernelIdentityResilienceToPositivity_Theorem :
  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by
  exact QFT_KernelIdentityConsolidationToPositivity_Theorem

