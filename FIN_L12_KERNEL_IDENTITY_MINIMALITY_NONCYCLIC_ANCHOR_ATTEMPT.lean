-- FIN Release 5.1: QW-2418 L12 kernel-identity-minimality noncyclic anchor attempt

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_KernelIdentityMinimalityToWellPosedness_Theorem :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  intro h
  exact RG_KernelIdentityClosureToWellPosedness_Theorem h
