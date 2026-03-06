# H30 Kernel-Invariant psi0 Anchor Audit

Status: `PASS_PARTIAL_KERNEL_INVARIANT_CANDIDATE_ONLY`
Date: `2026-03-06`

## Goal

Sprawdzic, czy `orientation_psi0 = f(phi, omega)` z dawnych modeli mozna uczciwie traktowac jako
wewnetrzny anchor orientacji wynikajacy z kernel invariants, czy nadal tylko jako parametr
preorientowanego kanalu anisotropowego.

## Inputs from prior studies

From `QW_1952` and `QW_1953`:
- `psi0 = mod(0.5*phi + 0.8*omega, 2*pi)`
- anisotropic term: `cos(2*(th-psi0)+phase)`
- old light/memory proxies: `retard_phase`, `observer_tau`, `observer_feedback_gain`

## Reduction

Two facts hold simultaneously.

1. `psi0` is deterministic from kernel-level quantities `phi` and `omega`.
   So it is stronger than a free narrative parameter.

2. In the old observer/light proxy models, `psi0` enters only as the orientation parameter of an
   already anisotropic channel.
   The models do not export:
   - a strict-core `theta_i` source,
   - a residual selector `2x2` block,
   - or a proof that `psi0` itself is the internal selector datum.

So the honest classification is:
- `psi0` is a kernel-invariant anchor candidate,
- but not yet a strict-core exported selector source.

## Result

The old proxy lane upgrades the status of `psi0` from a free heuristic angle to a deterministic
kernel-invariant candidate.

But this still falls short of the missing strict-core result.
There is no current export from `psi0` to actual `theta_1`, `theta_2`, and no proof that the
residual `O(2)` degeneracy is broken by `psi0` alone.

## Frontier

`H30_B1 := orientation_psi0 is a deterministic kernel-invariant anchor candidate derived from (phi,omega), but it is not yet exported or justified as the strict-core selector datum for theta_i and does not by itself discharge the residual O(2) degeneracy`

## Hard limits

- no `theorem-level PASS`
- no `full-closure PASS`
- no claim that `psi0` already equals the missing strict-core selector datum
- no claim that `psi0` exports actual `theta_1`, `theta_2`
- no claim that `QW-2191` is discharged
