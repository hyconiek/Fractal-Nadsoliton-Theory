-- QW-2442 strict attempt: single-foundation -> RG canonical export
set_option autoImplicit false
variable (NadsolitonSingleFoundation RG_WellPosedness_Target : Prop)

theorem RG_NadsolitonSingleFoundationToWellPosedness_EXPORT :
  NadsolitonSingleFoundation ->
  RG_WellPosedness_Target := by
  intro hFoundation
  exact RG_CanonicalAction_to_WellPosedness_EXPORT
