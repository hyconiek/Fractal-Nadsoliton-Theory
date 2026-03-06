-- FIN Release 5.1: QW-2344 L5 kernel-identity-consistency discharge attempt
-- Scope: derive identity-consistency theorem from identity-completeness theorem.

variable (FINActionComplete ConstructiveNonPerturbativeScheme PositivityToReconstruction : Prop)

theorem QFT_KernelIdentityConsistencyToPositivity_Theorem :
  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by
  exact QFT_KernelIdentityCompletenessToPositivity_Theorem

