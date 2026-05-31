# Phase-sign component-quotient dual-basis certificate

Status: `component-quotient-dual-basis-pairing-and-node-decomposition-certified-over-gf2`

## Result

This certificate exports the rows T_i and columns U_j as a Kronecker
dual basis and reconstructs the audited node vector from its coordinates.

## Summary

- pairing matrix identity: `True`
- rank pairing matrix: `12`
- expected coordinates: `True`
- active quotient coordinates only: `True`
- reconstructs node bits: `True`

## Hard limits

- K_strict_gate remains the current live/full operational kernel.
- No unqualified identity K_legacy_ont == K_strict_gate is claimed.
- No proof derives A(d), P(d), D(d), omega/phi, beta/eta, or the transport cocycle from strict nadsoliton dynamics.
- No beta_tors -> chi_11 theorem is claimed.
- No legacy physical-role transfer to K_strict_gate is claimed.
- No QW-2191 selector discharge is claimed.
- No ToE closure is claimed.
