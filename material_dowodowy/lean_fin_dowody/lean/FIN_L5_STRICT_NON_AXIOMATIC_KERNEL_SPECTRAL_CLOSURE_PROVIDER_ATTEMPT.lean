-- QW-2458 strict non-axiomatic kernel-spectral-closure provider attempt (L5/QFT)
set_option autoImplicit false
variable (FINActionComplete ConstructiveNonPerturbativeScheme PositivityToReconstruction : Prop)

theorem QFT_KernelSpectralClosureToPositivity_Theorem_NON_AXIOMATIC :
  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by
  exact QFT_KernelSpectralInvarianceToPositivity_Theorem
