-- FIN Release 5.1: QW-2308 L12 kernel-spectral-closure discharge attempt
-- Scope: derive kernel-spectral closure theorem from spectral-invariance theorem.

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_KernelSpectralClosureToWellPosedness_Theorem :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_KernelSpectralInvarianceToWellPosedness_Theorem

