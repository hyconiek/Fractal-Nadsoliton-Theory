import Mathlib

/-!
  FIN ST25 formal-core candidate.

  This source is intentionally restricted to algebraic identities that do not
  assert a physical interpretation.  It must not be called machine checked
  until it has been replayed with the pinned Lean/Mathlib toolchain recorded by
  the release.
-/

namespace FIN.ST25

theorem approximate_shadow_composition
    {R : Type*} [Ring R]
    (C₁ C₂ Φ Ψ Ξ R₁ R₂ : R)
    (h₁ : C₁ * Φ - Ψ * C₁ = R₁)
    (h₂ : C₂ * Ψ - Ξ * C₂ = R₂) :
    (C₂ * C₁) * Φ - Ξ * (C₂ * C₁) = C₂ * R₁ + R₂ * C₁ := by
  rw [← h₁, ← h₂]
  noncomm_ring

theorem three_vertex_dirichlet_identity
    {R : Type*} [CommRing R]
    (w₀₁ w₀₂ w₁₂ f₀ f₁ f₂ : R) :
    f₀ * (w₀₁ * (f₀ - f₁) + w₀₂ * (f₀ - f₂)) +
      f₁ * (w₀₁ * (f₁ - f₀) + w₁₂ * (f₁ - f₂)) +
      f₂ * (w₀₂ * (f₂ - f₀) + w₁₂ * (f₂ - f₁)) =
    w₀₁ * (f₀ - f₁)^2 + w₀₂ * (f₀ - f₂)^2 +
      w₁₂ * (f₁ - f₂)^2 := by
  ring

theorem thermal_scale_orbit
    {R : Type*} [Field R]
    (β E c : R) (hc : c ≠ 0) :
    (β / c) * (c * E) = β * E := by
  field_simp

theorem aut_z12_unit_orbit_contains_reverse :
    ((11 : ZMod 12) = -(1 : ZMod 12)) := by
  norm_num

end FIN.ST25
