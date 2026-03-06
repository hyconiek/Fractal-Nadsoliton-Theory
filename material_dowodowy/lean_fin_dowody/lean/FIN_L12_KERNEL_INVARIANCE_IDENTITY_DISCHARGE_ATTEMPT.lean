-- FIN Release 5.1: QW-2314 L12 kernel-invariance-identity discharge attempt
-- Scope: derive invariance-identity theorem from identity-minimality theorem.

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_KernelInvarianceIdentityToWellPosedness_Theorem :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_KernelIdentityMinimalityToWellPosedness_Theorem

