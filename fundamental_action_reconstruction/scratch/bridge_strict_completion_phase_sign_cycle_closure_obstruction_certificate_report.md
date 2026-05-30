# Phase-sign cycle-closure obstruction certificate

Status: `cycle-closure-parity-functional-certified-zero-closure-exact-odd-closure-obstructed`

## Result

The path sign cochain is tested against an artificial closing edge `11->0`.
The forced zero closing bit is exact on the 12-cycle; the odd closing-bit
perturbation is obstructed by the GF(2) cycle parity functional.

## Summary

- Rank(delta_cycle): `11`
- Nullity(delta_cycle): `1`
- H1 dimension: `1`
- Cycle rank: `1`
- Forced closing edge bit b11 xor b0: `0`
- Zero-closing cycle exact: `True`
- Zero-closing anchor recovers audited nodes: `True`
- Odd-closing cycle exact: `False`
- Odd-closing obstructed by cycle parity: `True`

## Hard limits

- K_strict_gate remains the current live/full operational kernel.
- No unqualified identity K_legacy_ont == K_strict_gate is claimed.
- No cyclic L5/L12 gate-generation strategy is introduced or recommended.
- No proof derives A(d), P(d), D(d), omega/phi, beta/eta, or the transport cocycle from strict nadsoliton dynamics.
- No beta_tors -> chi_11 theorem is claimed.
- No legacy physical-role transfer to K_strict_gate is claimed.
- No QW-2191 selector discharge is claimed.
- No ToE closure is claimed.
