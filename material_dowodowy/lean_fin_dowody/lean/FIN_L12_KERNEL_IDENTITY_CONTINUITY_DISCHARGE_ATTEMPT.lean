-- FIN Release 5.1: QW-2326 L12 kernel-identity-continuity discharge attempt
-- Scope: derive identity-continuity theorem from identity-coherence theorem.

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_KernelIdentityContinuityToWellPosedness_Theorem :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_KernelIdentityCoherenceToWellPosedness_Theorem

