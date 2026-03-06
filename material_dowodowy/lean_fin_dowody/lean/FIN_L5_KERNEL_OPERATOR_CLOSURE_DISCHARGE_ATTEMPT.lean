-- FIN Release 5.1: QW-2305 L5 kernel-operator-closure discharge attempt
-- Scope: derive kernel-operator theorem from kernel-spectral-closure theorem.

variable (FINActionComplete ConstructiveNonPerturbativeScheme PositivityToReconstruction : Prop)

theorem QFT_KernelOperatorClosureToPositivity_Theorem :
  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by
  exact QFT_KernelSpectralClosureToPositivity_Theorem

