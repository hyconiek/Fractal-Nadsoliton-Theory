-- FIN Release 5.1: QW-2371 L12 kernel-identity-universality discharge attempt
-- Scope: derive identity-universality theorem from identity-totality theorem.

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_KernelIdentityUniversalityToWellPosedness_Theorem :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_KernelIdentityTotalityToWellPosedness_Theorem

