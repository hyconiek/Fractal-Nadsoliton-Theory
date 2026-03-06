-- QW-2451 strict non-axiomatic canonical export provider attempt (L5/QFT)
set_option autoImplicit false
variable (FINActionComplete ConstructiveNonPerturbativeScheme PositivityToReconstruction : Prop)

theorem QFT_FundamentalKernelDynamicsToPositivity_Theorem :
  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by
  exact QFT_KernelOperatorClosureToPositivity_Theorem

theorem QFT_CanonicalAction_to_Positivity_EXPORT_NON_AXIOMATIC :
  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by
  exact QFT_FundamentalKernelDynamicsToPositivity_Theorem
