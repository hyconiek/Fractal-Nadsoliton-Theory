-- FIN Release 5.1: QW-2277 strict non-axiomatic residual provider attempt (RG)
-- File intentionally contains no explicit `axiom` declarations.

theorem RGGlobalWellPosednessAllScales_Derived :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_CanonicalAction_to_WellPosedness_EXPORT
