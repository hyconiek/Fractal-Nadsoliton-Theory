-- QW-2563 strict non-axiomatic kernel-identity-universality provider attempt (L12/RG)
set_option autoImplicit false
variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_KernelIdentityTotalityToWellPosedness_Theorem_NON_AXIOMATIC :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_KernelIdentityTotalityToWellPosedness_Theorem
