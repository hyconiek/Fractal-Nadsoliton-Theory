-- FIN Release 5.1: QW-2368 L5 kernel-identity-unification discharge attempt
-- Scope: derive identity-unification theorem from identity-universality theorem.

variable (FINActionComplete ConstructiveNonPerturbativeScheme PositivityToReconstruction : Prop)

theorem QFT_KernelIdentityUnificationToPositivity_Theorem :
  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by
  exact QFT_KernelIdentityUniversalityToPositivity_Theorem

