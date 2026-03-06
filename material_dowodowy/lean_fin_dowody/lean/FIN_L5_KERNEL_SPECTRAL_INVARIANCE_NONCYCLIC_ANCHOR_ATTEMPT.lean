-- FIN Release 5.1: QW-2410 L5 kernel-spectral-invariance noncyclic anchor attempt

variable (FINActionComplete ConstructiveNonPerturbativeScheme PositivityToReconstruction : Prop)

theorem QFT_KernelSpectralInvarianceToPositivity_Theorem :
  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by
  intro h
  exact QFT_KernelInvarianceIdentityToPositivity_Theorem h
