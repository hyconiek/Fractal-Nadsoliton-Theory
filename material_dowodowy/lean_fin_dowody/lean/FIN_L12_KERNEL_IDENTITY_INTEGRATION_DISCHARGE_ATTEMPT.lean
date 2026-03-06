-- FIN Release 5.1: QW-2365 L12 kernel-identity-integration discharge attempt
-- Scope: derive identity-integration theorem from identity-unification theorem.

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_KernelIdentityIntegrationToWellPosedness_Theorem :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_KernelIdentityUnificationToWellPosedness_Theorem

