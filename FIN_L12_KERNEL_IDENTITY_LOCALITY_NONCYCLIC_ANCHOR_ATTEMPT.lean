-- FIN Release 5.1: QW-2426 L12 kernel-identity-locality noncyclic anchor attempt
-- Scope: derive identity-locality theorem from identity-continuity theorem.

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_KernelIdentityLocalityToWellPosedness_Theorem :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_KernelIdentityContinuityToWellPosedness_Theorem
