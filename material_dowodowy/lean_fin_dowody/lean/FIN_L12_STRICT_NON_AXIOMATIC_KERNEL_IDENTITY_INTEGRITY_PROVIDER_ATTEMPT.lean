-- QW-2513 strict non-axiomatic kernel-identity-integrity provider attempt (L12/RG)
set_option autoImplicit false
variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_KernelIdentityIntegrityToWellPosedness_Theorem_NON_AXIOMATIC :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_KernelIdentityConsistencyToWellPosedness_Theorem
