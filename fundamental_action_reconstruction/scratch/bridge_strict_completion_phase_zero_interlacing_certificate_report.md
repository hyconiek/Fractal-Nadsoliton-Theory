# Strict completion phase-zero interlacing certificate

Status: `phase-transport-sign-flips-certified-by-legacy-strict-cosine-zero-interlacing`

## Result

This finite-domain audit explains the sign part of the phase/frequency
completion factor by exact cosine-zero interlacing on `[0,11]`.  It is a
certificate for consequences of the chosen phase parameters, not a strict
dynamical derivation of those parameters.

## Zero formula

`cos(omega*x+phi)=0 iff x=(pi/2+k*pi-phi)/omega`

## Interlacing summary

- Legacy zero positions: `[1.3333333333333335, 5.333333333333333, 9.333333333333332]`
- Strict zero positions: `[7.581676052731609]`
- Phase transport sign pattern: `[1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1]`
- Phase sign-flip edges: `['1->2', '5->6', '7->8', '9->10']`
- Zero parity matches all edges: `True`
- All integer nodes avoid phase zeros: `True`
- Min legacy zero endpoint separation: `3.333333333333e-01`
- Min strict zero endpoint separation: `4.183239472684e-01`

## Edge interlacing rows

| edge | legacy zeros inside | strict zeros inside | zero count | predicted flip | actual flip |
|---|---:|---:|---:|---:|---:|
| 0->1 | 0 | 0 | 0 | False | False |
| 1->2 | 1 | 0 | 1 | True | True |
| 2->3 | 0 | 0 | 0 | False | False |
| 3->4 | 0 | 0 | 0 | False | False |
| 4->5 | 0 | 0 | 0 | False | False |
| 5->6 | 1 | 0 | 1 | True | True |
| 6->7 | 0 | 0 | 0 | False | False |
| 7->8 | 0 | 1 | 1 | True | True |
| 8->9 | 0 | 0 | 0 | False | False |
| 9->10 | 1 | 0 | 1 | True | True |
| 10->11 | 0 | 0 | 0 | False | False |

## Guarded interpretation

The certificate explains the finite sign flips of the completion phase factor.
It does not derive strict omega/phi from nadsoliton dynamics, does not prove
`beta_tors -> chi_11`, and does not transfer legacy physical roles onto
`K_strict_gate`.

## Hard limits

- K_strict_gate remains the current live/full operational kernel.
- No unqualified identity K_legacy_ont == K_strict_gate is claimed.
- No proof derives strict omega/phi or phase transport from strict nadsoliton dynamics.
- No beta_tors -> chi_11 theorem is claimed.
- No legacy physical-role transfer to K_strict_gate is claimed.
- No QW-2191 selector discharge is claimed.
- No ToE closure is claimed.
