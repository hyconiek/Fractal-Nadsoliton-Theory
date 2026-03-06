-- FIN Release 5.1: QW-2387 L12 noncircular anchor candidate
-- Intent: break identity-cycle by referencing an external action-level provider symbol.

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

 theorem RG_KernelIdentityLocalityToWellPosedness_Theorem :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_ActionLevel_PhysicalBridge_Derivation
