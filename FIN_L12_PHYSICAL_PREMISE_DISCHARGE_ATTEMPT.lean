-- FIN Release 5.1: QW-2293 L12 physical-premise discharge execution attempt
-- Scope: axiom-token-free theorem-level attempt for RG physical bridge premise.

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_PhysicalBridgePremise :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_ActionLevel_PhysicalBridge_Derivation

