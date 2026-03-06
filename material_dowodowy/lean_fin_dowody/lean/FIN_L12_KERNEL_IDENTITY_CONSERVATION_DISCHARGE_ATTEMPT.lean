-- FIN Release 5.1: QW-2335 L12 kernel-identity-conservation discharge attempt
-- Scope: derive identity-conservation theorem from identity-compatibility theorem.

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_KernelIdentityConservationToWellPosedness_Theorem :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_KernelIdentityCompatibilityToWellPosedness_Theorem

