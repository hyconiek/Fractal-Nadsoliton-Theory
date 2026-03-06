-- FIN Release 5.1: QW-2329 L5 kernel-identity-coherence discharge attempt
-- Scope: derive identity-coherence theorem from identity-regularity theorem.

variable (FINActionComplete ConstructiveNonPerturbativeScheme PositivityToReconstruction : Prop)

theorem QFT_KernelIdentityCoherenceToPositivity_Theorem :
  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by
  exact QFT_KernelIdentityRegularityToPositivity_Theorem

