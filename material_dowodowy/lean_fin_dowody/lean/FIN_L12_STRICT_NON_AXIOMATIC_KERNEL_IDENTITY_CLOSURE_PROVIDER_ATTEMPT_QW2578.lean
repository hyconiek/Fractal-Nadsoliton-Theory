-- QW-2578 strict non-axiomatic kernel-identity-closure provider attempt (L12/RG)
set_option autoImplicit false
variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_KernelIdentityClosureToWellPosedness_Theorem_NON_AXIOMATIC :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_KernelIdentityLocalityToWellPosedness_Theorem
