-- FIN Release 5.1: QW-2335 L5 kernel-identity-conservation discharge attempt
-- Scope: derive identity-conservation theorem from identity-compatibility theorem.

variable (FINActionComplete ConstructiveNonPerturbativeScheme PositivityToReconstruction : Prop)

theorem QFT_KernelIdentityConservationToPositivity_Theorem :
  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by
  exact QFT_KernelIdentityCompatibilityToPositivity_Theorem

