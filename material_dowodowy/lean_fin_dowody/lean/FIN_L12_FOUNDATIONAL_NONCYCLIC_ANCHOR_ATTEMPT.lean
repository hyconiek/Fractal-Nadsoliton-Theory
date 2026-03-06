-- FIN Release 5.1: foundational noncyclic anchor candidate (L12)
-- Intent: move foundational frontier to fundamental-kernel-dynamics layer.

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

 theorem RG_FundamentalActionToWellPosedness_Derivation :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_FundamentalKernelDynamicsToWellPosedness_Theorem
