-- FIN Release 5.1: QW-2414 L5 kernel-invariance-identity noncyclic anchor attempt

variable (FINActionComplete ConstructiveNonPerturbativeScheme PositivityToReconstruction : Prop)

theorem QFT_KernelInvarianceIdentityToPositivity_Theorem :
  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by
  intro h
  exact QFT_KernelIdentityMinimalityToPositivity_Theorem h
