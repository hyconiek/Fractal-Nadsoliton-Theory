-- FIN Release 5.1: QW-2410 L12 kernel-spectral-invariance noncyclic anchor attempt

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_KernelSpectralInvarianceToWellPosedness_Theorem :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  intro h
  exact RG_KernelInvarianceIdentityToWellPosedness_Theorem h
