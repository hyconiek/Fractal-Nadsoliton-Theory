# Phase-zero carrier/edge incidence certificate

Status: `zero-carrier-to-edge-incidence-map-reconstructs-phase-z2-edge-bits`

## Result

This certificate writes rational zero-carrier crossings as a GF(2)
carrier-to-edge incidence matrix.  Multiplying by the all-one carrier
multiplicity vector recovers the audited phase edge bits.

## Summary

- Carrier order: `['legacy_z0', 'legacy_z1', 'strict_z0', 'legacy_z2']`
- Edge bits from incidence: `[0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0]`
- Flip edges from incidence: `['1->2', '5->6', '7->8', '9->10']`
- Column rank over GF(2): `4`
- All carriers strictly inside open edges: `True`
- All carrier columns single-edge standard basis: `True`
- All edges have at most one incident carrier: `True`
- Matches Z2 edge bits: `True`
- Matches GF(2) solution edge bits: `True`

## Edge incidence rows

| edge | incident carriers | edge bit | parity pass |
| --- | --- | ---: | :---: |
| 0->1 | [] | 0 | `True` |
| 1->2 | ['legacy_z0'] | 1 | `True` |
| 2->3 | [] | 0 | `True` |
| 3->4 | [] | 0 | `True` |
| 4->5 | [] | 0 | `True` |
| 5->6 | ['legacy_z1'] | 1 | `True` |
| 6->7 | [] | 0 | `True` |
| 7->8 | ['strict_z0'] | 1 | `True` |
| 8->9 | [] | 0 | `True` |
| 9->10 | ['legacy_z2'] | 1 | `True` |
| 10->11 | [] | 0 | `True` |

## Hard limits

- K_strict_gate remains the current live/full operational kernel.
- No unqualified identity K_legacy_ont == K_strict_gate is claimed.
- No proof derives A(d), P(d), D(d), omega/phi, beta/eta, or the transport cocycle from strict nadsoliton dynamics.
- No beta_tors -> chi_11 theorem is claimed.
- No legacy physical-role transfer to K_strict_gate is claimed.
- No QW-2191 selector discharge is claimed.
- No ToE closure is claimed.
