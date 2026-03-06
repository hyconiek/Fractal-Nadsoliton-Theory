-- FIN Release 5.1: QW-2350 L12 kernel-identity-saturation discharge attempt
-- Scope: derive identity-saturation theorem from identity-stability theorem.

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_KernelIdentitySaturationToWellPosedness_Theorem :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_KernelIdentityStabilityToWellPosedness_Theorem

