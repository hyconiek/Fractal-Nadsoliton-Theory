-- FIN Release 5.1: canonical export symbol bridge for RG (axiomatic layer only)
-- Scope boundary: this is NOT non-axiomatic closure; it is a symbol-completion bridge.

axiom FINActionComplete : Prop
axiom RGConstructiveMap : Prop
axiom RGGlobalWellPosednessAllScales : Prop

axiom RGGlobalWellPosednessAllScales_DerivedOrPending :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales

theorem RG_CanonicalAction_to_WellPosedness_EXPORT :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  intro h
  exact RGGlobalWellPosednessAllScales_DerivedOrPending h

