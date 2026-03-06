# RAPORT QW-2041: CANONICAL vs REFROZEN REPARAMETERIZATION AUDIT

- Data UTC: 2026-03-03T21:17:41.348395+00:00
- Verdict: **CANONICAL_REFROZEN_REPARAMETERIZATION_FAIL**
- Readiness: **CANONICAL_SEMANTIC_DRIFT_CONFIRMED_INTERNAL**

## Kernels
- Canonical TeX omega/phi/beta/eta: 0.785398 / 0.523599 / 0.010000 / 1.000000
- Refrozen QW-2039 omega/phi/beta/eta: 0.185750 / 0.162500 / 0.920000 / 1.800000

## Best Waveform Candidate (min RMSE)
- a/b/p/sign: 0.065618 / 1.101124 / 0.450000 / 1
- rmse: 0.057313
- corr: 0.673553
- r2: 0.434955
- affine_r2: 0.453674
- node errors median/q95/max rel: 2335.016927 / 5302.863651 / 5654.081817

## Best Gate Candidate (max hard flags)
- pass_count: 2/7
- a/b/p/sign: 0.065618 / 0.993258 / 0.587234 / 1
- rmse: 0.069517
- corr: 0.639703
- r2: 0.168686
- affine_r2: 0.409220
- node errors median/q95/max rel: 142.980586 / 233.790142 / 243.064664

## Gate Flags (best gate candidate)
- strict_corr_ge_0p95: False
- strict_r2_ge_0p90: False
- affine_r2_ge_0p95: False
- node_median_rel_le_0p10: False
- node_q95_rel_le_0p25: False
- physical_p_band_0p5_2p0: True
- physical_shift_absb_le_1: True

## Required Next Step
- DEFINE_EXPLICIT_BRIDGE_OPERATOR_BETWEEN_CANONICAL_AND_EFFECTIVE_SEMANTICS

## Artifacts
- JSON: `report_qw2041_canonical_refrozen_reparameterization_audit.json`
