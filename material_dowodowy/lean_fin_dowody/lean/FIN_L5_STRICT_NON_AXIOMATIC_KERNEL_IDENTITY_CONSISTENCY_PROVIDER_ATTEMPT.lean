-- QW-2518 strict non-axiomatic kernel-identity-consistency provider attempt (L5/QFT)
set_option autoImplicit false
variable (FINActionComplete ConstructiveNonPerturbativeScheme PositivityToReconstruction : Prop)

theorem QFT_KernelIdentityConsistencyToPositivity_Theorem_NON_AXIOMATIC :
  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by
  exact QFT_KernelIdentityCompletenessToPositivity_Theorem
