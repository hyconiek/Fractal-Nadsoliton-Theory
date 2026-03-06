-- FIN Release 5.1: QW-2317 L12 kernel-identity-minimality discharge attempt
-- Scope: derive identity-minimality theorem from identity-closure theorem.

variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)

theorem RG_KernelIdentityMinimalityToWellPosedness_Theorem :
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by
  exact RG_KernelIdentityClosureToWellPosedness_Theorem

