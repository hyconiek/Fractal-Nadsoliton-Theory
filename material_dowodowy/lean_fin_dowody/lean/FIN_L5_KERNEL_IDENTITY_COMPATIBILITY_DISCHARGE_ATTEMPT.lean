-- FIN Release 5.1: QW-2338 L5 kernel-identity-compatibility discharge attempt
-- Scope: derive identity-compatibility theorem from identity-integrity theorem.

variable (FINActionComplete ConstructiveNonPerturbativeScheme PositivityToReconstruction : Prop)

theorem QFT_KernelIdentityCompatibilityToPositivity_Theorem :
  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by
  exact QFT_KernelIdentityIntegrityToPositivity_Theorem

