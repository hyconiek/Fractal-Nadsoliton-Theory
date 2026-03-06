-- FIN Release 5.1: QW-2414 L12 kernel-invariance-identity noncyclic anchor attempt

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_KernelInvarianceIdentityToWellPosedness_Theorem :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  intro h
  exact RG_KernelIdentityMinimalityToWellPosedness_Theorem h
