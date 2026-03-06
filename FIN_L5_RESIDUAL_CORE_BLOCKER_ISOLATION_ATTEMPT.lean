-- FIN Release 5.1: QW-2282 QFT residual core-blocker isolation attempt
variable (FINActionComplete ConstructiveNonPerturbativeScheme PositivityToReconstruction : Prop)

theorem PositivityToReconstruction_Derived :
  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by
  exact QFT_CanonicalAction_to_Positivity_EXPORT
