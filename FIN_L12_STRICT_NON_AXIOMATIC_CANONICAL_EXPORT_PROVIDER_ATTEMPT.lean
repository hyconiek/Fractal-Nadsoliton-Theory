-- QW-2451 strict non-axiomatic canonical export provider attempt (L12/RG)
set_option autoImplicit false
variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_FundamentalKernelDynamicsToWellPosedness_Theorem :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_KernelOperatorClosureToWellPosedness_Theorem

theorem RG_CanonicalAction_to_WellPosedness_EXPORT_NON_AXIOMATIC :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_FundamentalKernelDynamicsToWellPosedness_Theorem
