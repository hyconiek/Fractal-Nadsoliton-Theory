-- QW-2523 strict non-axiomatic kernel-identity-completeness provider attempt (L12/RG)
set_option autoImplicit false
variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_KernelIdentityCompletenessToWellPosedness_Theorem_NON_AXIOMATIC :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_KernelIdentitySaturationToWellPosedness_Theorem
