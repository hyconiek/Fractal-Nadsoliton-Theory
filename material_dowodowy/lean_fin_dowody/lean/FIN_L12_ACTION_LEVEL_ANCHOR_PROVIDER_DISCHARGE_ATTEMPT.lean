-- FIN Release 5.1: QW-2390 L12 action-level anchor provider discharge attempt
-- Scope: build RG action-level provider from foundational derivation symbol.

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_ActionLevel_PhysicalBridge_Derivation :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_FundamentalActionToWellPosedness_Derivation

