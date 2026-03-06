-- FIN Release 5.1: QW-2422 L5 kernel-identity-closure noncyclic anchor attempt

variable (FINActionComplete ConstructiveNonPerturbativeScheme PositivityToReconstruction : Prop)

theorem QFT_KernelIdentityClosureToPositivity_Theorem :
  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by
  intro h
  exact QFT_KernelIdentityLocalityToPositivity_Theorem h
