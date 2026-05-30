# Strict completion phase-zero margin certificate

Status: `strict-phase-zero-edge-placement-has-positive-rational-parameter-margin`

## Result

This audit computes a rational robustness margin around the selected strict
phase parameters.  Within the certified symmetric parameter rectangle, the
strict k=0 zero remains in edge `7->8`, adjacent strict zeros remain outside
`[0,11]`, and the phase sign-flip pattern is preserved.

## Margin values

- k0 left margin: `22897/212000` = `1.080047169811e-01`
- k0 right margin: `17561/226000` = `7.770353982301e-02`
- k1 right-exclusion margin: `531381/212000` = `2.506514150943e+00`
- k-1 left-exclusion margin: `7349/4240` = `1.733254716981e+00`
- Limiting symmetric epsilon: `17561/2034000` = `8.633726647001e-03`
- Certified symmetric epsilon: `17561/4068000` = `4.316863323500e-03`

## Robustness summary

- All margins positive: `True`
- All worst-case inequalities hold at epsilon: `True`
- Active margin source: `k0_right_margin/9`
- Preserved flip edges: `['1->2', '5->6', '7->8', '9->10']`
- Preserved sign pattern: `[1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1]`

## Worst-case inequalities

| name | expression | slack | holds |
|---|---|---:|---:|
| k0_left_edge_upper_worst_case | `7*(omega+eps)+(phi+eps) < pi_lower/2` | 3168077/43120800 | True |
| k0_right_edge_lower_worst_case | `pi_upper/2 < 8*(omega-eps)+(phi-eps)` | 17561/452000 | True |
| k1_stays_right_of_d11_worst_case | `11*(omega+eps)+(phi+eps) < 3*pi_lower/2` | 176415227/71868000 | True |
| k_minus_1_stays_left_of_d0_worst_case | `-pi_lower/2 < phi-eps` | 372765917/215604000 | True |

## Guarded interpretation

This is only a local robustness certificate for the selected strict phase
parameters. It does not derive strict omega/phi from nadsoliton dynamics,
does not prove `beta_tors -> chi_11`, and does not transfer legacy physical
roles onto `K_strict_gate`.

## Hard limits

- K_strict_gate remains the current live/full operational kernel.
- No unqualified identity K_legacy_ont == K_strict_gate is claimed.
- No proof derives strict omega/phi or phase transport from strict nadsoliton dynamics.
- No beta_tors -> chi_11 theorem is claimed.
- No legacy physical-role transfer to K_strict_gate is claimed.
- No QW-2191 selector discharge is claimed.
- No ToE closure is claimed.
