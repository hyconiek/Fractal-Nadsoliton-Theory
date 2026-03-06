-- FIN Release 5.1: QW-2347 L12 kernel-identity-completeness discharge attempt
-- Scope: derive identity-completeness theorem from identity-saturation theorem.

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_KernelIdentityCompletenessToWellPosedness_Theorem :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_KernelIdentitySaturationToWellPosedness_Theorem

