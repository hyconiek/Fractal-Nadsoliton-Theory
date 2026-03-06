-- FIN Release 5.1: QW-2434 L12 kernel-identity-coherence noncyclic anchor attempt
-- Scope: derive identity-coherence theorem from identity-regularity theorem.

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_KernelIdentityCoherenceToWellPosedness_Theorem :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_KernelIdentityRegularityToWellPosedness_Theorem
