-- FIN Release 5.1: QW-2434 L5 kernel-identity-coherence noncyclic anchor attempt
-- Scope: derive identity-coherence theorem from identity-regularity theorem.

variable (FINActionComplete ConstructiveNonPerturbativeScheme PositivityToReconstruction : Prop)

theorem QFT_KernelIdentityCoherenceToPositivity_Theorem :
  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by
  exact QFT_KernelIdentityRegularityToPositivity_Theorem
