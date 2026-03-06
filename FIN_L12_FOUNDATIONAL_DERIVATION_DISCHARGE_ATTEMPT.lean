-- FIN Release 5.1: QW-2299 L12 foundational derivation discharge attempt
-- Scope: derive foundational RG bridge from fundamental kernel dynamics theorem.

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_FundamentalActionToWellPosedness_Derivation :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_FundamentalKernelDynamicsToWellPosedness_Theorem

