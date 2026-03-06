-- FIN Release 5.1: QW-2356 L5 kernel-identity-robustness discharge attempt
-- Scope: derive identity-robustness theorem from identity-resilience theorem.

variable (FINActionComplete ConstructiveNonPerturbativeScheme PositivityToReconstruction : Prop)

theorem QFT_KernelIdentityRobustnessToPositivity_Theorem :
  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by
  exact QFT_KernelIdentityResilienceToPositivity_Theorem

