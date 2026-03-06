-- FIN Release 5.1: QW-2377 L12 kernel-identity-finality discharge attempt
-- Scope: derive identity-finality theorem from identity-closure theorem.

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_KernelIdentityFinalityToWellPosedness_Theorem :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_KernelIdentityClosureToWellPosedness_Theorem

