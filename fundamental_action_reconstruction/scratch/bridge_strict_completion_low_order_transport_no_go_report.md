# Strict completion low-order transport no-go audit

Status: `low-order-autonomous-transport-readings-fail-for-completion-cocycle`

## Result

This finite-domain audit tests whether the completion transport cocycle can
already be explained by very simple autonomous laws.  All tested low-order
readings fail on `d=0..11`; this sharpens but does not close the open
`strict_transport_derivation` blocker.

## Transport input summary

- Edge sign pattern: `[1, -1, 1, 1, 1, -1, 1, -1, 1, -1, 1]`
- Negative edges: `['1->2', '5->6', '7->8', '9->10']`
- Edge ratio spread: `9.928523012792e+00`

## No-go summary

- Positive damping-only transport fails: `True`
- Constant first-order multiplier fails: `True`
- Constant second-order recurrence fails: `True`
- Affine log-envelope fails: `True`
- Short-period edge-sign law fails: `True`
- Minimum failed-model L2 residual: `5.666440484247e-02`

## Model residuals

- First-order best `r`: `2.930754805735e-01`, L2 residual `6.376000586350e-01`
- Second-order best `(a,b)`: `(-1.590684005367e-01, -7.646138093067e-02)`, L2 residual `5.666440484247e-02`
- Affine log best `(m,c)`: `(-5.099555766966e-01, -1.218092366137e+00)`, L2 residual `3.634601695668e+00`

## Guarded interpretation

The result is a finite no-go for the listed low-order/autonomous readings only.
It is not a full no-go theorem against strict nadsoliton dynamics and does not
derive any selector, chi_11 source, or ToE closure.

## Hard limits

- K_strict_gate remains the current live/full operational kernel.
- No unqualified identity K_legacy_ont == K_strict_gate is claimed.
- No full no-go theorem against every strict transport derivation is claimed.
- No proof derives the completion transport from strict nadsoliton dynamics.
- No beta_tors -> chi_11 theorem is claimed.
- No legacy physical-role transfer to K_strict_gate is claimed.
- No QW-2191 selector discharge is claimed.
- No ToE closure is claimed.
