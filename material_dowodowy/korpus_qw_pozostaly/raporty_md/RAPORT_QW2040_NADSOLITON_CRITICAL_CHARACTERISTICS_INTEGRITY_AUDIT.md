# RAPORT QW-2040: NADSOLITON CRITICAL CHARACTERISTICS INTEGRITY AUDIT

- Data UTC: 2026-03-03T20:57:57.927274+00:00

## Refrozen 2039 Kernel
- omega/phi/beta/eta: 0.185750 / 0.162500 / 0.920000 / 1.800000

## Canonical-Formula Drift vs TeX (beta_tors semantics)
- alpha_EM^-1 rel_err: 9.991e-01
- n_s rel_err: 1.857e+00
- G_ratio_20 rel_err: 1.887e+14
- k_tau rel_err: 6.849e+00
- zero_spacing rel_err: 3.228e+00
- first_zero rel_err: 4.686e+00

## TeX Canonical Verdict: **TEX_CANONICAL_CHARACTERISTICS_NOT_PRESERVED**
## Current Refrozen Branch Verdict: **CURRENT_REFROZEN_BRANCH_CHARACTERISTICS_OPERATIONALLY_PRESERVED**

## Interpretation
- Canonical TeX semantics are not preserved numerically in refrozen branch, but current refrozen branch remains internally consistent and gate-valid.

## Artifacts
- JSON: `report_qw2040_nadsoliton_critical_characteristics_integrity_audit.json`
