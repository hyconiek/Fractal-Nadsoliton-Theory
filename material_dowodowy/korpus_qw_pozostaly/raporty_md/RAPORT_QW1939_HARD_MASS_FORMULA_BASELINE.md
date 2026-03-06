# RAPORT QW-1939: HARD MASS FORMULA BASELINE

- Data UTC: 2026-03-03T06:19:30.080367+00:00
- Kernel: omega=0.373414, phi=-1.310234, beta=0.615938, eta=2.80
- Verdict: **HARD_MASS_FORMULA_BASELINE_FAIL**

## Exact Formula
- m(Q) = m_top * 4^(-(gamma * Q / 4))
- No Delta correction, no fitting.

## Variants
- canonical_gamma_1p52: gamma=1.520000, mean/max rel%=8.242/18.904
- kernel_derived_gamma_1to2: gamma=2.395104, mean/max rel%=78.792/99.925
- kernel_derived_gamma_1to4: gamma=2.349948, mean/max rel%=78.152/99.890

## Primary Hard Baseline
- variant: kernel_derived_gamma_1to4
- hard_pass: False

## Required Next Step
- MASS_SECTOR_REMAINS_OPEN_UNDER_EXACT_FORMULA

## Artifacts
- JSON: `report_qw1939_hard_mass_formula_baseline.json`
