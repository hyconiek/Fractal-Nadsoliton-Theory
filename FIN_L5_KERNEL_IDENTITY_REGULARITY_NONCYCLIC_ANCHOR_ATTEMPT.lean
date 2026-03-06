-- FIN Release 5.1: QW-2438 L5 kernel-identity-regularity noncyclic anchor attempt
-- Scope: derive identity-regularity theorem from identity-conservation theorem.

variable (FINActionComplete ConstructiveNonPerturbativeScheme PositivityToReconstruction : Prop)

theorem QFT_KernelIdentityRegularityToPositivity_Theorem :
  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by
  exact QFT_KernelIdentityConservationToPositivity_Theorem
