-- FIN Release 5.1: QW-2380 L5 kernel-identity-closure discharge attempt
-- Scope: derive identity-closure theorem from identity-locality theorem.

variable (FINActionComplete ConstructiveNonPerturbativeScheme PositivityToReconstruction : Prop)

theorem QFT_KernelIdentityClosureToPositivity_Theorem :
  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by
  exact QFT_KernelIdentityLocalityToPositivity_Theorem

