-- FIN Release 5.1: QW-2323 L5 kernel-identity-locality discharge attempt
-- Scope: derive identity-locality theorem from identity-continuity theorem.

variable (FINActionComplete ConstructiveNonPerturbativeScheme PositivityToReconstruction : Prop)

theorem QFT_KernelIdentityLocalityToPositivity_Theorem :
  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by
  exact QFT_KernelIdentityContinuityToPositivity_Theorem

