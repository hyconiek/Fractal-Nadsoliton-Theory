-- FIN Release 5.1: QW-2374 L12 kernel-identity-totality discharge attempt
-- Scope: derive identity-totality theorem from identity-finality theorem.

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_KernelIdentityTotalityToWellPosedness_Theorem :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_KernelIdentityFinalityToWellPosedness_Theorem

