-- QW-2442 strict attempt: single-foundation -> QFT canonical export
set_option autoImplicit false
variable (NadsolitonSingleFoundation QFT_Positivity_Target : Prop)

theorem QFT_NadsolitonSingleFoundationToPositivity_EXPORT :
  NadsolitonSingleFoundation ->
  QFT_Positivity_Target := by
  intro hFoundation
  exact QFT_CanonicalAction_to_Positivity_EXPORT
