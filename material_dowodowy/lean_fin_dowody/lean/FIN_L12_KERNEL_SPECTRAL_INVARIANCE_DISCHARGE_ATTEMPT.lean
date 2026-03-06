-- FIN Release 5.1: QW-2311 L12 kernel-spectral-invariance discharge attempt
-- Scope: derive spectral-invariance theorem from invariance-identity theorem.

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_KernelSpectralInvarianceToWellPosedness_Theorem :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_KernelInvarianceIdentityToWellPosedness_Theorem

