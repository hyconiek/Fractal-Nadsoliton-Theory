-- FIN Release 5.1: QW-2380 L12 kernel-identity-closure discharge attempt
-- Scope: derive identity-closure theorem from identity-locality theorem.

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_KernelIdentityClosureToWellPosedness_Theorem :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_KernelIdentityLocalityToWellPosedness_Theorem

