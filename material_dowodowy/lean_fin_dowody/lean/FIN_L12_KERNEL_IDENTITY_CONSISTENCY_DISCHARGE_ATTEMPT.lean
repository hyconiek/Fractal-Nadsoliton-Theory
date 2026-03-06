-- FIN Release 5.1: QW-2344 L12 kernel-identity-consistency discharge attempt
-- Scope: derive identity-consistency theorem from identity-completeness theorem.

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_KernelIdentityConsistencyToWellPosedness_Theorem :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_KernelIdentityCompletenessToWellPosedness_Theorem

