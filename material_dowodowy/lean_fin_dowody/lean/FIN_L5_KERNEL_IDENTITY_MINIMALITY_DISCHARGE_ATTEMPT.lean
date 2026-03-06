-- FIN Release 5.1: QW-2317 L5 kernel-identity-minimality discharge attempt
-- Scope: derive identity-minimality theorem from identity-closure theorem.

variable (FINActionComplete ConstructiveNonPerturbativeScheme PositivityToReconstruction : Prop)

theorem QFT_KernelIdentityMinimalityToPositivity_Theorem :
  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by
  exact QFT_KernelIdentityClosureToPositivity_Theorem

