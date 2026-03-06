-- FIN Release 5.1: QW-2341 L12 kernel-identity-integrity discharge attempt
-- Scope: derive identity-integrity theorem from identity-consistency theorem.

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_KernelIdentityIntegrityToWellPosedness_Theorem :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_KernelIdentityConsistencyToWellPosedness_Theorem

