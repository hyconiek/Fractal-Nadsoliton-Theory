# V6 Psi0+Viscosity Distinct Effect Audit

Status: `PASS_PARTIAL_SPECTRAL_SPLIT_ONLY`
Date: `2026-03-06`

## Goal

Sprawdzic, czy lane `psi0 + anizotropowa viscosity` daje nowy pair-level efekt,
ktorego nie ma juz w czystym lane `psi0`, czy tylko przeskalowuje ten sam
zaimportowany kierunek.

## Inputs

- `H31`: coordinate-level embedding `u_psi0_pair1 = cos(psi0)c_1 + sin(psi0)s_1`
- `V3`: admissible anisotropic pair-level viscosity class
- `V4`: coupled operator `K_visc_aniso^(1)(psi0) = R(psi0) diag(nu_parallel, nu_perp) R(-psi0)`
- `V5`: anti-overclaim boundary certificate for `psi0 + viscosity`

## Comparison

Pure `psi0` lane supplies only a direction candidate on `V_1`:

`u_psi0_pair1 = cos(psi0)c_1 + sin(psi0)s_1`

Coupled `psi0 + viscosity` lane adds one thing only:

- spectral weighting / splitting between:
  - `u_parallel(psi0)`
  - `u_perp(psi0)`

So the coupled lane may generate:
- amplitude asymmetry,
- damping asymmetry,
- response-time asymmetry,

but only relative to an already imported anchor direction `psi0`.

## Result

The honest result is:
- `psi0 + viscosity` does produce a pair-level effect not present in the bare coordinate embedding alone,
- but the new effect is only a spectral/response split around the already given anchor,
- not a new orientation datum.

So the coupled lane is:
- stronger than bare `psi0` as a response model,
- not stronger than bare `psi0` as a selector-source model.

## Frontier

`V6_B1 := psi0_plus_viscosity contributes a genuine spectral/response split on V_1 beyond the bare psi0 coordinate embedding, but it still imports rather than generates the anchor and therefore does not provide a new selector source`

## Hard limits

- no `theorem-level PASS`
- no `full-closure PASS`
- no claim that viscosity creates new orientation information
- no claim that the coupled lane supersedes the primary `psi0` anchor lane
- no claim that `QW-2191` is discharged
