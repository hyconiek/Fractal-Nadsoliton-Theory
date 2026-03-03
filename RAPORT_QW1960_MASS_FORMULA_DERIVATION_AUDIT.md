# RAPORT QW-1960: MASS FORMULA DERIVATION AUDIT

- Data UTC: 2026-03-03T07:53:03.615090+00:00
- Verdict: **DERIVATION_CONTAINS_MATERIAL_ERRORS_AND_CIRCULAR_STEPS**

## Canonical Formula Check (gamma=1.52)
- mean rel error: 8.242%
- max rel error: 18.904%
- implied gamma mean/std from fixed Q table: 1.5268 / 0.0344

## Arithmetic Check of Documented Gamma Origin
- expression tested: gamma = 2.26 / 1.77
- computed gamma: 1.276836
- delta vs 1.52: -0.243164
- if this gamma is used in mass law: mean/max rel error = 213.736% / 726.482%

## Circularity Check
- mean |Q_inferred - Q_model| (non-top): 0.185
- infer Q from masses (gamma=1.52) -> round Q -> reconstruct masses: mean/max rel error = 8.242% / 18.904%
- interpretation: low error here is not independent confirmation.

## Structural Degeneracy Check
- tau and charm share Q=9 in current mapping.
- experimental tau/charm ratio: 1.3991
- tau/charm gap: 39.913%
- formula prediction at Q=9: 1510.083 MeV (single value for both)

## Frozen-Kernel Compatibility
- gamma from kernel (1->4): 2.349948
- delta vs 1.52: +0.829948

## Flags
- arithmetic_inconsistency_gamma_origin: True
- q_mapping_circularity_risk: True
- tau_charm_same_q_non_degeneracy: True
- frozen_kernel_gamma_incompatibility: True

## Required Next Step
- RE-DERIVE_GAMMA_AND_Q_FROM_INDEPENDENT_TOPOLOGY_WITHOUT_MASS_INVERSION

## Artifacts
- JSON: `report_qw1960_mass_formula_derivation_audit.json`
