-- FIN Release 5.1: QW-2326 L5 kernel-identity-continuity discharge attempt
-- Scope: derive identity-continuity theorem from identity-coherence theorem.

variable (FINActionComplete ConstructiveNonPerturbativeScheme PositivityToReconstruction : Prop)

theorem QFT_KernelIdentityContinuityToPositivity_Theorem :
  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by
  exact QFT_KernelIdentityCoherenceToPositivity_Theorem

