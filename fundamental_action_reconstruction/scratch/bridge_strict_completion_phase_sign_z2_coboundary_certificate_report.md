# Strict completion phase-sign Z2 coboundary certificate

Status: `phase-sign-pattern-certified-as-z2-coboundary-on-finite-path-domain`

## Result

This audit packages the certified phase sign pattern as a finite Z2
cochain/coboundary on the path graph `0--1--...--11`.  It verifies
prefix reconstruction and all interval parity laws without introducing a
new cosine fit, selector premise, or strict dynamical derivation.

## Z2 summary

- Node bit pattern: `[0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0]`
- Edge bit pattern: `[0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0]`
- Flip support size: `4`
- Derived flip edges: `['1->2', '5->6', '7->8', '9->10']`
- Interval count: `66`
- Prefix reconstruction passes: `True`
- All interval coboundary rows pass: `True`

## Edge bit rows

| edge | left bit | right bit | edge bit | flip |
|---|---:|---:|---:|---:|
| 0->1 | 0 | 0 | 0 | False |
| 1->2 | 0 | 1 | 1 | True |
| 2->3 | 1 | 1 | 0 | False |
| 3->4 | 1 | 1 | 0 | False |
| 4->5 | 1 | 1 | 0 | False |
| 5->6 | 1 | 0 | 1 | True |
| 6->7 | 0 | 0 | 0 | False |
| 7->8 | 0 | 1 | 1 | True |
| 8->9 | 1 | 1 | 0 | False |
| 9->10 | 1 | 0 | 1 | True |
| 10->11 | 0 | 0 | 0 | False |

## Proof certificate

- `node_step`: Convert the certified phase signs into Z2 node bits b(d).
- `edge_step`: Compute edge bits as the coboundary e(d,d+1)=b(d) xor b(d+1).
- `prefix_step`: Starting from b(0), prefix xor of edge bits reconstructs every node bit and phase sign.
- `interval_step`: All 66 nontrivial intervals satisfy sum(edge bits) mod 2 = endpoint node-bit xor.
- `nonduplication`: This is a Z2 cochain/coboundary algebra certificate, not another zero-location, cell-sign, damping, or real-valued cocycle audit.
- `theoretical_limit`: The certificate proves finite sign-transport algebra on the selected domain; it does not derive omega/phi or transport from strict nadsoliton dynamics and does not discharge a selector obstruction.

## Hard limits

- K_strict_gate remains the current live/full operational kernel.
- No unqualified identity K_legacy_ont == K_strict_gate is claimed.
- No proof derives strict omega/phi or phase transport from strict nadsoliton dynamics.
- No beta_tors -> chi_11 theorem is claimed.
- No legacy physical-role transfer to K_strict_gate is claimed.
- No QW-2191 selector discharge is claimed.
- No ToE closure is claimed.
