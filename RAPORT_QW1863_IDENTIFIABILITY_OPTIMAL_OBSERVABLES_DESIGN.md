# RAPORT QW-1863: IDENTIFIABILITY OPTIMAL OBSERVABLES DESIGN

- Data UTC: 2026-03-03T00:04:01.023353+00:00
- Verdict: **IDENTIFIABILITY_DESIGN_PARTIAL**

## Target Tuple
- omega=0.785398
- phi=0.523599
- beta=0.010000

## Best Baseline Feature Set
- features: phase_increment, envelope_decay, zero_cross_offset
- score: 0.5242
- cond_q90: 5.699e+00
- beta-omega coupling (median): 0.080

## Signed Intervention Result
- score: 0.5247
- cond_q90: 5.688e+00
- beta-omega coupling (median): 0.079

## Gains
- score ratio (intervention/base): 1.001
- cond_q90 reduction factor: 1.002

## Top Baseline Combos
- k=3, score=0.5242, cond_q90=5.699e+00, coupling=0.080 | phase_increment, envelope_decay, zero_cross_offset
- k=3, score=0.5035, cond_q90=6.281e+00, coupling=0.095 | phase_increment, zero_cross_offset, envelope_kurtosis_proxy
- k=5, score=0.4868, cond_q90=8.179e+00, coupling=0.069 | phase_increment, zero_cross_offset, signed_asymmetry, torsion_cross_term, phase_curvature
- k=4, score=0.4857, cond_q90=9.354e+00, coupling=0.043 | phase_increment, zero_cross_offset, signed_asymmetry, torsion_cross_term
- k=4, score=0.4836, cond_q90=8.202e+00, coupling=0.074 | phase_increment, zero_cross_offset, signed_asymmetry, envelope_kurtosis_proxy
- k=3, score=0.4802, cond_q90=9.059e+00, coupling=0.060 | phase_increment, zero_cross_offset, torsion_cross_term
- k=5, score=0.4704, cond_q90=8.271e+00, coupling=0.098 | phase_increment, zero_cross_offset, signed_asymmetry, phase_curvature, envelope_kurtosis_proxy
- k=5, score=0.4702, cond_q90=8.687e+00, coupling=0.088 | phase_increment, envelope_decay, zero_cross_offset, torsion_cross_term, phase_curvature
- k=5, score=0.4638, cond_q90=8.796e+00, coupling=0.098 | phase_increment, zero_cross_offset, torsion_cross_term, phase_curvature, envelope_kurtosis_proxy
- k=5, score=0.4632, cond_q90=8.524e+00, coupling=0.106 | phase_increment, envelope_decay, zero_cross_offset, phase_curvature, envelope_kurtosis_proxy
- k=4, score=0.4596, cond_q90=1.069e+01, coupling=0.067 | phase_increment, envelope_decay, zero_cross_offset, signed_asymmetry
- k=4, score=0.4594, cond_q90=8.837e+00, coupling=0.106 | phase_increment, envelope_decay, zero_cross_offset, phase_curvature

## Recommended Next Studies
- QW-1865: Synthetic recovery benchmark on best observable set -> Verify unbiased recovery of beta, omega, phi under realistic noise.
- QW-1866: Paired signed-intervention execution protocol -> Implement baseline vs sign-flipped branch measurements to break omega-beta coupling.

## Artifacts
- JSON: `report_qw1863_identifiability_optimal_observables_design.json`
