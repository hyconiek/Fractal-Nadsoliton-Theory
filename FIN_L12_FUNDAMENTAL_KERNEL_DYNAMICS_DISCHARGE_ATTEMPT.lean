-- FIN Release 5.1: QW-2302 L12 fundamental-kernel-dynamics discharge attempt
-- Scope: derive kernel dynamics theorem from kernel-operator-closure theorem.

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_FundamentalKernelDynamicsToWellPosedness_Theorem :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_KernelOperatorClosureToWellPosedness_Theorem

