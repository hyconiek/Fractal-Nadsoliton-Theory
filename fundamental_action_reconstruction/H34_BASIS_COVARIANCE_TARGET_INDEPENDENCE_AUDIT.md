# H34 Basis Covariance / Target Independence Audit

Status: `PASS_PARTIAL_NO_STRICT_COVARIANCE_ARGUMENT`
Date: `2026-03-06`

## Goal

Sprawdzic, czy obecny strict core zawiera jakikolwiek argument typu:

- `basis covariance`
- `target independence`

ktory pozwalalby traktowac redukcje `psi0` nie jako wybor lokalnego chartu `pair1`,
ale jako fizycznie stabilna redukcje niezalezna od chartu.

## Inputs

- `H30`: `psi0` jako kernel-invariant anchor candidate.
- `H31`: legalny embedding `psi0 -> pair1`.
- `H33`: `pair1` jest tylko deterministic local chart, nie selector-relevant target.
- `QW-2190`: mode scaffold i deterministic pairs.
- `QW-2191`: residualna degeneracja `O(2)` pozostaje nierozladowana.

## Audit

Current strict-core support gives:

- deterministic local charts,
- deterministic mode pairs,
- legal local coordinate embeddings,
- local formulas on each `span{c_i,s_i}`.

Current strict-core support does **not** give:

- a theorem that embeddings related by local `O(2)` basis change are selector-equivalent,
- a theorem that reducing `psi0` on `pair1` is physically equivalent to reducing it on another admissible pair,
- a target-independence statement that would lift the embedding beyond chart dependence,
- a covariance export from the residual selector lane.

So the current honest classification is:

- local chart machinery: present,
- strict covariance/target-independence argument: absent.

## Result

The repository currently supports only a local-chart embedding story for `psi0`.

It does not yet support the stronger statement that the `psi0` reduction is
chart-independent or basis-covariant in a way sufficient for selector closure.

## Frontier

`H34_B1 := strict core contains local chart embeddings for psi0 but no basis-covariance or target-independence argument elevating those embeddings beyond chart dependence`

## Hard limits

- no `theorem-level PASS`
- no `full-closure PASS`
- no claim that `psi0 -> pair1` is chart-independent
- no claim that target choice is irrelevant
- no claim that `QW-2191` is discharged
