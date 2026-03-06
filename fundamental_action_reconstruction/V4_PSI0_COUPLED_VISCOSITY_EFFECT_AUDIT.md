# V4 Psi0-Coupled Viscosity Effect Audit

Status: `PASS_PARTIAL_COMPETING_EXTENSION_HYPOTHESIS_ONLY`
Date: `2026-03-06`

## Goal

Sprawdzic, czy sprzezenie:

- kernel-invariant anchor candidate `psi0`,
- z anizotropowa para-level viscosity class na `V_1 = span{c_1,s_1}`

produkuje nowy, niezalezny selector mechanism,

czy tylko modyfikuje juz zaimportowany anchor.

## Coupled operator

Take the admissible anisotropic viscosity class from `V3`:

`K_visc_aniso^(1)(psi0) = R(psi0) diag(nu_parallel, nu_perp) R(-psi0)`

with:

- `nu_parallel >= 0`
- `nu_perp >= 0`
- possibly `nu_parallel != nu_perp`

and with `psi0` imported as the anchor candidate from `H30/H31`.

## Pair-level consequence

On `V_1`, this coupled operator has:

- preferred eigenvector: `u_parallel(psi0) = cos(psi0) c_1 + sin(psi0) s_1`
- orthogonal eigenvector: `u_perp(psi0) = -sin(psi0) c_1 + cos(psi0) s_1`

So the operator can:

- amplify, damp, or split response relative to the anchor direction,
- but it does **not** generate `psi0` by itself.

## Interpretation

This means the coupled effect is:

- `anchor-amplifying`
- or `anchor-refining`

not:

- `anchor-generating`

Therefore the coupled lane is meaningful as a secondary closure aid,

but it does not replace the primary need for an anchor source.

## Result

The honest current result is:

- `psi0 + anisotropic viscosity` yields a nontrivial pair-level operator,
- but the nontriviality is conditional on an already supplied anchor,
- so this lane does not supersede `psi0` and does not by itself solve selector closure.

## Frontier

`V4_B1 := coupling anisotropic viscosity to psi0 yields a nontrivial pair-level anchor-amplifying operator on V_1, but it does not generate the anchor itself and therefore cannot by itself close the selector problem`

## Hard limits

- no `theorem-level PASS`
- no `full-closure PASS`
- no claim that viscosity generates `psi0`
- no claim that the coupled lane replaces the primary `psi0` anchor lane
- no claim that `QW-2191` is discharged
