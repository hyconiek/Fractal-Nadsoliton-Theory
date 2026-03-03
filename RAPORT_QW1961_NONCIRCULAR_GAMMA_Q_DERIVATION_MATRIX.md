# RAPORT QW-1961: NONCIRCULAR GAMMA/Q DERIVATION MATRIX

- Data UTC: 2026-03-03T08:00:28.258898+00:00
- Verdict: **NONCIRCULAR_DERIVATION_HAS_STRICT_PASS_CANDIDATE**
- variants total: 16
- noncircular pass count: 1

## Key Inputs
- D_f = 4 ln 2 = 2.772589
- n_grav = 2.2600
- gamma_kernel(d1->d4) = 2.349948
- delta_info (from QW-1958 lambda/mu) = 0.177027

## Best Overall Variant
- q/gamma/split: legacy_fibonacci / canonical_frozen_1p52_reference=1.520000 / info_split_from_qw1958
- mean/max err: 8.880% / 20.044%
- tau/charm ratio pred/exp/error: 1.2050 / 1.3991 / 13.872%
- all_pass: True

## Best Noncircular Variant
- q/gamma/split: legacy_fibonacci / derived_force_energy_2n_over_3=1.506667 / info_split_from_qw1958
- mean/max err: 12.051% / 34.013%
- tau/charm ratio pred/exp/error: 1.2031 / 1.3991 / 14.013%
- all_pass: True

## Required Next Step
- PROMOTE_BEST_NONCIRCULAR_BRANCH_TO_UNIFIED_TRIAD_TEST

## Artifacts
- JSON: `report_qw1961_noncircular_gamma_q_derivation_matrix.json`
