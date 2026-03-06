# H31 psi0-to-pair1 Reduction Audit

Status: `PASS_PARTIAL_COORDINATE_EMBEDDING_ONLY`
Date: `2026-03-06`

## Goal

Sprawdzic, czy `orientation_psi0` mozna uczciwie zredukowac do aktualnej pary
`pair1 = (c_1, s_1)` bez przemycania selekcji orientacji.

## Inputs

- `H30`: `psi0 = mod(0.5*phi + 0.8*omega, 2*pi)` jest deterministycznym kandydatem anchoru z kernel invariants.
- `C3`: `pair1 = (c_1, s_1)` jest aktualna deterministyczna para modowa.
- `C34`: dla dowolnego kata lokalny reprezentant ma postac
  `u_1(theta_1) = cos(theta_1)c_1 + sin(theta_1)s_1`.

## Reduction

The formally available embedding is:

`u_psi0_pair1 := cos(psi0) c_1 + sin(psi0) s_1`

This is algebraically legal inside:

`V_1 = span{c_1, s_1}`

But three gaps remain:

1. no strict-core export says that `pair1` is the selector-relevant reduction target for `psi0`;
2. no strict-core rule says that `psi0` should be interpreted as `theta_1` on `pair1`;
3. no basis-covariance argument shows that this embedding is more than a coordinate choice.

## Result

There is a coordinate-level candidate reduction:

`psi0 -> u_psi0_pair1`

but not yet a provenance-valid strict-core reduction:

`psi0 -> theta_1 on pair1`

So the current honest classification is:
- coordinate embedding candidate: present,
- strict-core selector reduction: absent.

## Frontier

`H31_B1 := psi0 admits a formal coordinate-level embedding into pair1 via u_psi0_pair1 = cos(psi0)c_1 + sin(psi0)s_1, but there is no strict-core justification that this embedding is the selector-relevant reduction rather than a coordinate choice`

## Hard limits

- no `theorem-level PASS`
- no `full-closure PASS`
- no claim that `psi0 = theta_1`
- no claim that `pair1` is already proven to be the selector-relevant reduction target
- no claim that `QW-2191` is discharged
