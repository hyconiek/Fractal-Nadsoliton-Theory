-- FIN Release 5.1: QW-2281 RG residual core-blocker isolation attempt
variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RGGlobalWellPosednessAllScales_Derived :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_CanonicalAction_to_WellPosedness_EXPORT
