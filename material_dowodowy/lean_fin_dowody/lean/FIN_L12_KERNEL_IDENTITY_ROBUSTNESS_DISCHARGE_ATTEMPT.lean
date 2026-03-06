-- FIN Release 5.1: QW-2356 L12 kernel-identity-robustness discharge attempt
-- Scope: derive identity-robustness theorem from identity-resilience theorem.

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_KernelIdentityRobustnessToWellPosedness_Theorem :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_KernelIdentityResilienceToWellPosedness_Theorem

