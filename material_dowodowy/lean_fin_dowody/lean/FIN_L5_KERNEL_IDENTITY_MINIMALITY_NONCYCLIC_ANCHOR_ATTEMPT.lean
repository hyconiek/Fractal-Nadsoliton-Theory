-- FIN Release 5.1: QW-2418 L5 kernel-identity-minimality noncyclic anchor attempt

variable (FINActionComplete ConstructiveNonPerturbativeScheme PositivityToReconstruction : Prop)

theorem QFT_KernelIdentityMinimalityToPositivity_Theorem :
  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by
  intro h
  exact QFT_KernelIdentityClosureToPositivity_Theorem h
