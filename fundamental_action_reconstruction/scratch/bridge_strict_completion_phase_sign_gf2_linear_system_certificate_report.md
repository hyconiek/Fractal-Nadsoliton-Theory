# Phase-sign GF(2) linear-system certificate

Status: `gf2-prefix-system-full-rank-with-unique-four-edge-phase-flip-solution`

## Result

This certificate writes finite prefix reconstruction as a GF(2) linear
system `L e = b_tail xor b0`.  The prefix matrix is full-rank, has an
explicit first-difference inverse, and has exactly the audited four-flip
solution.

## Summary

- Rank: `11`
- Nullity: `0`
- Determinant mod 2: `1`
- Unique solution: `True`
- Solution Hamming weight: `4`
- Flip edges: `['1->2', '5->6', '7->8', '9->10']`
- All equations pass: `True`
- Inherits edge-support minimality: `True`

## Equation rows

| node d | rhs | evaluated | passes |
| ---: | ---: | ---: | :---: |
| 1 | 0 | 0 | `True` |
| 2 | 1 | 1 | `True` |
| 3 | 1 | 1 | `True` |
| 4 | 1 | 1 | `True` |
| 5 | 1 | 1 | `True` |
| 6 | 0 | 0 | `True` |
| 7 | 0 | 0 | `True` |
| 8 | 1 | 1 | `True` |
| 9 | 1 | 1 | `True` |
| 10 | 0 | 0 | `True` |
| 11 | 0 | 0 | `True` |

## Hard limits

- K_strict_gate remains the current live/full operational kernel.
- No unqualified identity K_legacy_ont == K_strict_gate is claimed.
- No proof derives A(d), P(d), D(d), omega/phi, beta/eta, or the transport cocycle from strict nadsoliton dynamics.
- No beta_tors -> chi_11 theorem is claimed.
- No legacy physical-role transfer to K_strict_gate is claimed.
- No QW-2191 selector discharge is claimed.
- No ToE closure is claimed.
