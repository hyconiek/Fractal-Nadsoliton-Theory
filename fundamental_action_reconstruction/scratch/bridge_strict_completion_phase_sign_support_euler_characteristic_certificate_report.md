# Phase-sign support Euler-characteristic certificate

Status: `support-induced-graph-euler-characteristic-and-boundary-count-certified`

## Result

This certificate computes the Euler characteristic of the induced
1-node support graph and checks the corresponding path-boundary count.

## Summary

- Support nodes: `[2, 3, 4, 5, 8, 9]`
- Internal edges: `['2->3', '3->4', '4->5', '8->9']`
- Components: `[[2, 3, 4, 5], [8, 9]]`
- Euler V-E: `6 - 4 = 2`
- Component count: `2`
- Boundary weight: `4`
- Boundary edges equal flip support: `True`

## Hard limits

- K_strict_gate remains the current live/full operational kernel.
- No unqualified identity K_legacy_ont == K_strict_gate is claimed.
- No proof derives A(d), P(d), D(d), omega/phi, beta/eta, or the transport cocycle from strict nadsoliton dynamics.
- No beta_tors -> chi_11 theorem is claimed.
- No legacy physical-role transfer to K_strict_gate is claimed.
- No QW-2191 selector discharge is claimed.
- No ToE closure is claimed.
