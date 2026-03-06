-- FIN Release 5.1: QW-2341 L5 kernel-identity-integrity discharge attempt
-- Scope: derive identity-integrity theorem from identity-consistency theorem.

variable (FINActionComplete ConstructiveNonPerturbativeScheme PositivityToReconstruction : Prop)

theorem QFT_KernelIdentityIntegrityToPositivity_Theorem :
  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by
  exact QFT_KernelIdentityConsistencyToPositivity_Theorem

