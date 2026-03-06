-- FIN Release 5.1: QW-2377 L5 kernel-identity-finality discharge attempt
-- Scope: derive identity-finality theorem from identity-closure theorem.

variable (FINActionComplete ConstructiveNonPerturbativeScheme PositivityToReconstruction : Prop)

theorem QFT_KernelIdentityFinalityToPositivity_Theorem :
  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by
  exact QFT_KernelIdentityClosureToPositivity_Theorem

