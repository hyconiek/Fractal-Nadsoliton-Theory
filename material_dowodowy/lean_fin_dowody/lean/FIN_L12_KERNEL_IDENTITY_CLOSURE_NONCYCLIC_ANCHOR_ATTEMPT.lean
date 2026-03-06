-- FIN Release 5.1: QW-2422 L12 kernel-identity-closure noncyclic anchor attempt

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_KernelIdentityClosureToWellPosedness_Theorem :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  intro h
  exact RG_KernelIdentityLocalityToWellPosedness_Theorem h
