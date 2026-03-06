-- FIN Release 5.1: QW-2314 L5 kernel-invariance-identity discharge attempt
-- Scope: derive invariance-identity theorem from identity-minimality theorem.

variable (FINActionComplete ConstructiveNonPerturbativeScheme PositivityToReconstruction : Prop)

theorem QFT_KernelInvarianceIdentityToPositivity_Theorem :
  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by
  exact QFT_KernelIdentityMinimalityToPositivity_Theorem

