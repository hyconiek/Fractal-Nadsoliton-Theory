-- FIN Release 5.1: QW-2368 L12 kernel-identity-unification discharge attempt
-- Scope: derive identity-unification theorem from identity-universality theorem.

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_KernelIdentityUnificationToWellPosedness_Theorem :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_KernelIdentityUniversalityToWellPosedness_Theorem

