-- FIN Release 5.1: QW-2296 L12 action-level provider discharge attempt
-- Scope: build RG action-level provider from foundational derivation symbol.

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_ActionLevel_PhysicalBridge_Derivation :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_FundamentalActionToWellPosedness_Derivation

