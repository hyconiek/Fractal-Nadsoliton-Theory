-- FIN Release 5.1: QW-2353 L5 kernel-identity-stability discharge attempt
-- Scope: derive identity-stability theorem from identity-robustness theorem.

variable (FINActionComplete ConstructiveNonPerturbativeScheme PositivityToReconstruction : Prop)

theorem QFT_KernelIdentityStabilityToPositivity_Theorem :
  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by
  exact QFT_KernelIdentityRobustnessToPositivity_Theorem

