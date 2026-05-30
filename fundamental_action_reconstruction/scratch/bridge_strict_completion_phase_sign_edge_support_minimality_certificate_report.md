# Phase-sign edge-support minimality certificate

Status: `unique-four-edge-z2-support-is-minimal-for-the-audited-phase-sign-pattern`

## Result

This certificate treats the phase-sign Z2 data as a finite path-graph
problem.  It exhaustively enumerates all `2^11` edge-bit assignments and
checks which assignments reconstruct the audited node-bit pattern from the
left anchor.

## Summary

- Total edge assignments checked: `2048`
- Matching assignment count: `1`
- Unique matching assignment: `True`
- Target support size: `4`
- Flip edges: `['1->2', '5->6', '7->8', '9->10']`
- Lower-support assignments checked: `232`
- All lower-support assignments fail: `True`

## Lower-support exhaustion

| support size | checked | all fail |
| ---: | ---: | :---: |
| 0 | 1 | `True` |
| 1 | 11 | `True` |
| 2 | 55 | `True` |
| 3 | 165 | `True` |

## Hard limits

- K_strict_gate remains the current live/full operational kernel.
- No unqualified identity K_legacy_ont == K_strict_gate is claimed.
- No proof derives A(d), P(d), D(d), omega/phi, beta/eta, or the transport cocycle from strict nadsoliton dynamics.
- No beta_tors -> chi_11 theorem is claimed.
- No legacy physical-role transfer to K_strict_gate is claimed.
- No QW-2191 selector discharge is claimed.
- No ToE closure is claimed.
