-- FIN Release 5.1: QW-2305 L12 kernel-operator-closure discharge attempt
-- Scope: derive kernel-operator theorem from kernel-spectral-closure theorem.

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_KernelOperatorClosureToWellPosedness_Theorem :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_KernelSpectralClosureToWellPosedness_Theorem

