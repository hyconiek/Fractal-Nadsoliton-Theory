-- FIN Release 5.1: QW-2338 L12 kernel-identity-compatibility discharge attempt
-- Scope: derive identity-compatibility theorem from identity-integrity theorem.

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_KernelIdentityCompatibilityToWellPosedness_Theorem :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_KernelIdentityIntegrityToWellPosedness_Theorem

