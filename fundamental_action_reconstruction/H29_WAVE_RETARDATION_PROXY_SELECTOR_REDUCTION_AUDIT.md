# H29 Wave/Retardation Proxy Selector Reduction Audit

Status: `PASS_PARTIAL_PREORIENTED_PROXY_REDUCTION_ONLY`
Date: `2026-03-06`

## Goal

Sprawdzic, czy stare proxy falowo-opoznieniowe:
- `retard_phase`
- `observer_tau`
- `observer_feedback_gain`

same moga wygenerowac selektor orientacji, czy tylko moduluja kanal, ktory ma juz wstepnie ustawiona orientacje.

## Inputs from prior studies

From `QW_1952` / `QW_1953` / `QW_1956`:
- anisotropic term: `cos(2*(th-psi0)+phase)`
- phase proxy: `retard_phase`
- memory proxy: `observer_tau` or `tau_h/tau_l`
- gain proxy: `observer_feedback_gain`
- orientation anchor: `orientation_psi0`

## Reduction

Expand the anisotropic term:

`cos(2*(th-psi0)+phase) = cos(2*(th-psi0))*cos(phase) - sin(2*(th-psi0))*sin(phase)`

This shows:
- `retard_phase` changes weights of already present angular harmonics,
- `observer_tau` changes time/memory scale,
- `observer_feedback_gain` changes amplitude,
- but `orientation_psi0` still supplies the directional anchor.

## Result

The old wave/retardation proxies do **not** by themselves generate an internal strict-core orientation anchor.
They only modulate a pre-oriented anisotropic channel.

So the most honest conclusion is:
- existing light/memory proxies are relevant,
- but they do not yet solve the `theta_i` selector problem,
- and they do not collapse `QW-2191` on their own.

## Frontier

`H29_B1 := existing wave/retardation/memory proxies only modulate a pre-oriented anisotropic channel and do not by themselves generate an internal strict-core orientation anchor for theta_i`

## Hard limits

- no `theorem-level PASS`
- no `full-closure PASS`
- no claim that `orientation_psi0` is already strict-core internalized
- no claim that old light/observer calculations already solve the selector degeneracy
