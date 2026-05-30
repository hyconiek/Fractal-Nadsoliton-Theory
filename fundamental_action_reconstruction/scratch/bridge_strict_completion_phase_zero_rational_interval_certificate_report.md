# Strict completion phase-zero rational interval certificate

Status: `phase-zero-sign-flip-edges-certified-by-rational-pi-intervals`

## Result

This audit strengthens the phase-zero sign-flip certificate by using rational
interval arithmetic instead of relying on floating zero placement for the
strict phase zero.

## Rational inputs

- Pi interval: `333/106 < pi < 355/113`
- Strict omega: `743/4000`
- Strict phi: `13/80`
- Exact legacy zeros: `['4/3', '16/3', '28/3']`

## Interval summary

- Strict k=-1 left of domain: `True`
- Strict k=0 inside edge 7->8: `True`
- Strict k=1 right of audit domain: `True`
- Strict zero count in (0,11): `1`
- Legacy zero count in (0,11): `3`
- Flip edges from rational intervals: `['1->2', '5->6', '7->8', '9->10']`
- Sign pattern from rational intervals: `[1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1]`
- Matches float phase-zero report flip edges: `True`
- Matches float phase-zero report sign pattern: `True`

## Strict zero interval rows

| k | x lower | x upper | placement |
|---:|---:|---:|---|
| -1 | -783450/83959 | -367450/39379 | left of 0 |
| 0 | 298550/39379 | 636550/83959 | inside 7->8 |
| 1 | 964550/39379 | 2056550/83959 | right of 11 |

## Guarded interpretation

This is only a rational-interval certificate for phase-zero placement. It does
not derive strict omega/phi from nadsoliton dynamics, does not prove
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
