# P709 Current Strict Release‑7 OS Residual‑Sign Gauge‑Irrelevance Audit Probe (No False‑PASS)

Status: `P709_CURRENT_STRICT_RELEASE_7_OS_RESIDUAL_SIGN_GAUGE_IRRELEVANCE_AUDIT_PROBE_SPEC_NO_FALSE_PASS`  
As of: `2026-03-17`

## Goal

`T173` keeps the strict post‑`T172` frontier honest:

- kernel‑alone/global `QW-2191` discharge remains unclaimed, and
- any directed/sign‑sensitive physical orientation datum in strict core remains unclaimed.

However, Release‑7 also exports an **Operational ToE closure** in **projective** scope (`N701`) and multiple downstream
OS observables (`P694`, `P696`, `F704`, …) that are *intended* to be sign‑gauge‑safe.

This probe audits the minimal concrete claim needed to continue honestly without a false “need to break `QW-2191` first”:

> For the Release‑7 OS observables actually used downstream (quadratic proxy spectra and invariant mass observable),
> the residual per‑pair sign choice `u_m -> -u_m` does not change the exported numeric results (within tolerance).

This is an **audit** only; it does not upgrade any directed sign into strict core.

## Inputs

1. Projective closure object:
   - `generated/selector_closure_global_c_v1_projective_strict_v1.json` (`F672`)
2. Diagonal/local Hessian value instantiation:
   - `generated/psi_hessian_diagonal_local_strict_derived_value_instantiated_v1.json` (`F459`)
3. Release‑7 OS outputs (baseline for comparison):
   - `generated/p694_*_summary.json` (`P694`)
   - `generated/p696_*_summary.json` (`P696`)
   - `generated/f704_*_summary.json` (`F704`)

## Acceptance (no false pass)

For each sign pattern `(s_1,..,s_5)∈{±1}^5` applied to the exported `pair1..pair5` rays:

1. `P694` invariants are unchanged:
   - `m2_proxy(pair_m)` (Rayleigh quotient on symmetrized `H_psi`) matches baseline within tolerance.
2. `P696` invariants are unchanged:
   - diagonal channel proxy spectrum `diag(B^T H_psi B)` matches baseline within tolerance,
   - mixing scalars reported as Frobenius‑norm ratios (`offdiag_to_diag_ratio`, `offblock_max_fro`) match baseline.
3. `F704` invariants are unchanged:
   - `min/max m2_proxy` from the eigenvalue spectrum of `H_psi` match baseline within tolerance.

## Hard limits

This probe does **not**:

- discharge kernel‑alone/global `QW-2191`,
- discharge `T173`,
- claim any directed/sign‑sensitive physical orientation datum in strict core,
- claim Standard Model identification,
- claim physical‑unit calibration (proxy→GeV),
- claim “actual emergent observer closure”.

