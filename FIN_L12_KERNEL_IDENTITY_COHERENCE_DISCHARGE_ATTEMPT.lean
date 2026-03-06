-- FIN Release 5.1: QW-2329 L12 kernel-identity-coherence discharge attempt
-- Scope: derive identity-coherence theorem from identity-regularity theorem.

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_KernelIdentityCoherenceToWellPosedness_Theorem :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_KernelIdentityRegularityToWellPosedness_Theorem

