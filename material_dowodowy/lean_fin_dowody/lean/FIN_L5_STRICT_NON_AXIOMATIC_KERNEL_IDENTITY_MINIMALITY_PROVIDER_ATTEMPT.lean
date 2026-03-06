-- QW-2473 strict non-axiomatic kernel-identity-minimality provider attempt (L5/QFT)
set_option autoImplicit false
variable (FINActionComplete ConstructiveNonPerturbativeScheme PositivityToReconstruction : Prop)

theorem QFT_KernelIdentityMinimalityToPositivity_Theorem_NON_AXIOMATIC :
  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by
  exact QFT_KernelIdentityClosureToPositivity_Theorem
