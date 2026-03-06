-- FIN Release 5.1: QW-2308 L5 kernel-spectral-closure discharge attempt
-- Scope: derive kernel-spectral closure theorem from spectral-invariance theorem.

variable (FINActionComplete ConstructiveNonPerturbativeScheme PositivityToReconstruction : Prop)

theorem QFT_KernelSpectralClosureToPositivity_Theorem :
  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by
  exact QFT_KernelSpectralInvarianceToPositivity_Theorem

