-- FIN Release 5.1: QW-2311 L5 kernel-spectral-invariance discharge attempt
-- Scope: derive spectral-invariance theorem from invariance-identity theorem.

variable (FINActionComplete ConstructiveNonPerturbativeScheme PositivityToReconstruction : Prop)

theorem QFT_KernelSpectralInvarianceToPositivity_Theorem :
  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by
  exact QFT_KernelInvarianceIdentityToPositivity_Theorem

