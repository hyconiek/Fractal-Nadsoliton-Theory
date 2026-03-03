# RAPORT QW-1729: NADSOLITON -> KERNEL CHARACTERISTICS MAP

- Data UTC: 2026-03-02T16:34:59.080358+00:00
- Kernel closure score: 0.522
- Werdykt: **KERNEL_CHARACTERISTICS_NOT_CLOSED**

## Mapowanie charakterystyk nadsolitona na kernel
- amplitude_alpha_geo: 4-bit informational capacity / info-geometry identity -> alpha_geo | status=postulated_consistency_supported | score=0.65
- oscillation_frequency: octave resonance periodicity / 8-octave cycle -> omega = pi/4 | status=geometric_ansatz | score=0.50
- phase_offset: hexagonal lattice/vacuum symmetry -> phi = pi/6 | status=geometric_ansatz | score=0.50
- damping_scale: inter-layer torsion damping -> beta_tors = 0.01 | status=calibrated_then_frozen | score=0.45
- hyperbolic_denominator: fractal topology + topological tunneling path summation -> 1/(1 + beta_tors*d) | status=mechanistic_reduction | score=0.55
- combined_effective_form: joint effect of info capacity, resonance, torsion, topology -> K(d)=alpha_geo*cos(omega*d+phi)/(1+beta_tors*d) | status=effective_model | score=0.48

## Nierozwiazane punkty
- Microscopic derivation of beta_tors (not only calibration).
- First-principles derivation of (omega, phi) from Nadsoliton dynamics.
- Formal proof of hyperbolic denominator from path integral on fractal manifold.
- Flavor-sector derivation (CKM/PMNS) from the same kernel with shared parameters.

## Artefakty
- JSON szczegolowy: `report_qw1729_nadsoliton_kernel_characteristics_map.json`
