-- FIN Release 5.1: QW-2359 L12 kernel-identity-resilience discharge attempt
-- Scope: derive identity-resilience theorem from identity-consolidation theorem.

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_KernelIdentityResilienceToWellPosedness_Theorem :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_KernelIdentityConsolidationToWellPosedness_Theorem

