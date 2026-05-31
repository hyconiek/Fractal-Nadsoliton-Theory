# Phase-sign node-support interval-boundary certificate

Status: `node-support-maximal-interval-boundary-recovers-edge-flip-cochain`

## Result

This certificate decomposes the audited 1-node support into maximal path
intervals and verifies that their GF(2) boundary is exactly the edge-flip
cochain.

## Summary

- Maximal intervals: `[{'start': 2, 'end': 5}, {'start': 8, 'end': 9}]`
- Interval-boundary rank: `2`
- Interval-boundary nullity: `0`
- Node bits recovered from interval union: `True`
- Edge bits recovered from interval boundaries: `True`
- Boundary weight formula: `2*2 - 0 = 4`
- Flip edges: `['1->2', '5->6', '7->8', '9->10']`

## Hard limits

- K_strict_gate remains the current live/full operational kernel.
- No unqualified identity K_legacy_ont == K_strict_gate is claimed.
- No proof derives A(d), P(d), D(d), omega/phi, beta/eta, or the transport cocycle from strict nadsoliton dynamics.
- No beta_tors -> chi_11 theorem is claimed.
- No legacy physical-role transfer to K_strict_gate is claimed.
- No QW-2191 selector discharge is claimed.
- No ToE closure is claimed.
