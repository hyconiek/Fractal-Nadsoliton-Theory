# Phase-zero carrier-prefix node-matrix certificate

Status: `carrier-prefix-node-matrix-recovers-cell-sign-and-z2-node-bits`

## Result

This certificate composes the GF(2) prefix matrix with the carrier/edge
incidence matrix, producing a carrier-prefix node matrix `C=L*M` whose
rows list the zero-carriers left of each integer node.

## Summary

- Carrier order: `['legacy_z0', 'legacy_z1', 'strict_z0', 'legacy_z2']`
- Node bits from carrier-prefix matrix: `[0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0]`
- Phase signs from carrier-prefix matrix: `[1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1]`
- Carrier-prefix rank over GF(2): `4`
- All node rows match cell-sign left carriers: `True`
- All node bits match cell-sign: `True`
- Edge differences recover carrier/edge incidence: `True`
- Matches Z2 node bits: `True`

## Node prefix rows

| d | carriers left by matrix | node bit | sign | matches cell-sign |
| ---: | --- | ---: | ---: | :---: |
| 0 | [] | 0 | 1 | `True` |
| 1 | [] | 0 | 1 | `True` |
| 2 | ['legacy_z0'] | 1 | -1 | `True` |
| 3 | ['legacy_z0'] | 1 | -1 | `True` |
| 4 | ['legacy_z0'] | 1 | -1 | `True` |
| 5 | ['legacy_z0'] | 1 | -1 | `True` |
| 6 | ['legacy_z0', 'legacy_z1'] | 0 | 1 | `True` |
| 7 | ['legacy_z0', 'legacy_z1'] | 0 | 1 | `True` |
| 8 | ['legacy_z0', 'legacy_z1', 'strict_z0'] | 1 | -1 | `True` |
| 9 | ['legacy_z0', 'legacy_z1', 'strict_z0'] | 1 | -1 | `True` |
| 10 | ['legacy_z0', 'legacy_z1', 'strict_z0', 'legacy_z2'] | 0 | 1 | `True` |
| 11 | ['legacy_z0', 'legacy_z1', 'strict_z0', 'legacy_z2'] | 0 | 1 | `True` |

## Hard limits

- K_strict_gate remains the current live/full operational kernel.
- No unqualified identity K_legacy_ont == K_strict_gate is claimed.
- No proof derives A(d), P(d), D(d), omega/phi, beta/eta, or the transport cocycle from strict nadsoliton dynamics.
- No beta_tors -> chi_11 theorem is claimed.
- No legacy physical-role transfer to K_strict_gate is claimed.
- No QW-2191 selector discharge is claimed.
- No ToE closure is claimed.
