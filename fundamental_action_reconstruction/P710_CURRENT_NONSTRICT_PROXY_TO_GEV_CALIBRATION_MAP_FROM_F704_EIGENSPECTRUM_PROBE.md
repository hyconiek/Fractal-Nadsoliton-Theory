# P710 Current Nonstrict Proxy$\to$GeV Calibration Map From `F704` Eigen‑Spectrum Probe (No False‑PASS)

Status: `P710_CURRENT_NONSTRICT_PROXY_TO_GEV_CALIBRATION_MAP_FROM_F704_EIGENSPECTRUM_PROBE_SPEC_NO_FALSE_PASS`  
As of: `2026-03-17`

## Goal

Release‑7 strict scope exports only **dimensionless** quadratic proxy spectra (`P694`, `P696`) and the basis‑invariant
eigen‑spectrum mass observable (`F704`), with meaning discipline (`N703`) and with **no strict physical‑unit map**.

To make the theory operationally testable against external experimental numbers, any proxy$\to$GeV mapping must be:

1. explicitly **non‑strict**,
2. split into an explicit **external dataset** and an explicit **policy** object,
3. frozen and versioned to prevent silent degrees‑of‑freedom smuggling,
4. and, unless pass criteria are explicitly defined, must default to **diagnostic-only** (no pass claim).

This probe exports one such non‑strict calibration map for the basis‑invariant `F704` eigenvalue spectrum using a single
global scale parameter (no sectors, no kinetic refit, no nonlinear remap).

## Inputs

1. Proxy observable (strict, dimensionless):
   - `generated/mass_observable_diagonal_local_strict_derived_v1.json` (`F704`)
2. External dataset (nonstrict):
   - `external_data/sm_mass_targets_v1.json`
3. Calibration policy (nonstrict):
   - `external_data/proxy_to_gev_calibration_policy_v1.json`

## Output

1. `scale_m2_per_proxy_unit` (GeV$^2$ per proxy unit), and `scale_m_per_sqrt_proxy_unit` (GeV per $\sqrt{\text{proxy}}$),
2. a calibrated proxy spectrum in GeV (no identification / no assignment),
3. hard limits explicit: no Standard Model match claim and no strict-unit claim.

