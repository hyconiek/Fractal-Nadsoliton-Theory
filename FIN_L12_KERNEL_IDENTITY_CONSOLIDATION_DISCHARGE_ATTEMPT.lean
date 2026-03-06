-- FIN Release 5.1: QW-2362 L12 kernel-identity-consolidation discharge attempt
-- Scope: derive identity-consolidation theorem from identity-integration theorem.

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_KernelIdentityConsolidationToWellPosedness_Theorem :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_KernelIdentityIntegrationToWellPosedness_Theorem

