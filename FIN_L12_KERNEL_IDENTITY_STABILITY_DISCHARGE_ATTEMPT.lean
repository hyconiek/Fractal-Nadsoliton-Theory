-- FIN Release 5.1: QW-2353 L12 kernel-identity-stability discharge attempt
-- Scope: derive identity-stability theorem from identity-robustness theorem.

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_KernelIdentityStabilityToWellPosedness_Theorem :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_KernelIdentityRobustnessToWellPosedness_Theorem

