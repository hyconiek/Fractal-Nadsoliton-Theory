-- FIN Release 5.1: QW-2332 L12 kernel-identity-regularity discharge attempt
-- Scope: derive identity-regularity theorem from identity-conservation theorem.

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_KernelIdentityRegularityToWellPosedness_Theorem :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_KernelIdentityConservationToWellPosedness_Theorem

