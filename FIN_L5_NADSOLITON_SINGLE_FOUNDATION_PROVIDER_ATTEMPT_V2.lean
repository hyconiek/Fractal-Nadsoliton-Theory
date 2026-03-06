-- QW-2445 strict attempt V2: single-foundation -> QFT canonical export
set_option autoImplicit false
variable (NadsolitonSingleFoundation QFT_Positivity_Target : Prop)

theorem QFT_NadsolitonSingleFoundationToPositivity_EXPORT_V2 :
  NadsolitonSingleFoundation ->
  QFT_Positivity_Target := by
  intro hFoundation
  exact QFT_CanonicalAction_to_Positivity_EXPORT
