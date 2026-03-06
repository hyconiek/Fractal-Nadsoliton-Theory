-- FIN Release 5.1: QW-2438 L12 kernel-identity-regularity noncyclic anchor attempt
-- Scope: derive identity-regularity theorem from identity-conservation theorem.

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_KernelIdentityRegularityToWellPosedness_Theorem :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_KernelIdentityConservationToWellPosedness_Theorem
