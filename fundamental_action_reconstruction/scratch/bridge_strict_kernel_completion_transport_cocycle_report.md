# Strict Kernel Completion Transport Cocycle Certificate

Status: `unique-finite-node-transport-and-adjacent-edge-cocycle-certified-for-completion-ansatz`

## Result

This finite-domain audit turns the pointwise completion quotient into a
unique node transport `T(d)` and adjacent-edge cocycle `R(d->d+1)` on
`d=0..11`.  It is a transport witness induced by the completion ansatz,
not a derivation of the ansatz from strict nadsoliton dynamics.

## Definitions

- `node_transport`: T(d)=K_strict_gate(d)/K_legacy_ont(d)=A(d)*P(d)*D(d)
- `edge_cocycle`: R(d->d+1)=T(d+1)/T(d)
- `additive_log_split`: Delta log|T| = Delta log|P| + Delta log D; Delta log A = 0 on edges
- `sign_cocycle`: sign R records phase-sign flips forced by the strict-vs-legacy cosine transport

## Cocycle summary

- Max node factorization residual: `1.110223024625e-16`
- Max reconstruction residual from edges: `8.673617379884e-19`
- Max edge log-split residual: `4.440892098501e-16`
- Interval count: `66`
- Max interval multiplicative residual: `8.881784197001e-16`
- Max interval additive log residual: `1.776356839400e-15`
- Transport sign pattern: `[1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1]`
- Phase sign-flip edges: `['1->2', '5->6', '7->8', '9->10']`
- `T(0)`: `4.109835665709e-01`
- `T(11)`: `3.238803935929e-03`
- `T(11)/T(0)`: `7.880616645945e-03`

## Edge cocycle rows

| edge | R | sign | log_abs_R | delta_log_abs_phase | delta_log_damping | split_residual |
|---|---:|---:|---:|---:|---:|---:|
| 0->1 | 1.609532879848e+00 | 1 | 4.759440001582e-01 | 1.159140849865e+00 | -6.831968497068e-01 | -3.330669073875e-16 |
| 1->2 | -2.136092891097e-01 | -1 | -1.543606684500e+00 | -7.464916614878e-01 | -7.971150230124e-01 | 0.000000000000e+00 |
| 2->3 | 2.488546082716e-01 | 1 | -1.390886455579e+00 | -7.936184631421e-01 | -5.972679924365e-01 | -2.220446049250e-16 |
| 3->4 | 5.792878584490e-01 | 1 | -5.459557601031e-01 | -8.818171952351e-02 | -4.577740405796e-01 | 2.220446049250e-16 |
| 4->5 | 1.733414391600e+00 | 1 | 5.500931002102e-01 | 9.166573471957e-01 | -3.665642469855e-01 | -1.110223024625e-16 |
| 5->6 | -2.397414553656e-01 | -1 | -1.428194205622e+00 | -1.124234748239e+00 | -3.039594573826e-01 | 0.000000000000e+00 |
| 6->7 | 1.488078243372e-01 | 1 | -1.905099575054e+00 | -1.646324381647e+00 | -2.587751934070e-01 | 2.220446049250e-16 |
| 7->8 | -6.412520101475e-01 | -1 | -4.443327477691e-01 | -2.195448393129e-01 | -2.247879084563e-01 | -1.665334536938e-16 |
| 8->9 | 9.205568343771e+00 | 1 | 2.219808555720e+00 | 2.418171759204e+00 | -1.983632034845e-01 | 4.440892098501e-16 |
| 9->10 | -7.229546690204e-01 | -1 | -3.244087572393e-01 | -1.471450360551e-01 | -1.772637211842e-01 | 1.110223024625e-16 |
| 10->11 | 6.024741049754e-01 | 1 | -5.067105938348e-01 | -3.466647197574e-01 | -1.600458740774e-01 | 0.000000000000e+00 |

## Guarded interpretation

The cocycle is the unique finite transport object forced after choosing the
completion ansatz. It does not prove that strict nadsoliton dynamics generates
that transport, does not prove `beta_tors -> chi_11`, and does not transfer
legacy physical roles onto `K_strict_gate`.

## Hard limits

- K_strict_gate remains the current live/full operational kernel.
- No unqualified identity K_legacy_ont == K_strict_gate is claimed.
- No proof derives the transport cocycle from strict nadsoliton dynamics.
- No beta_tors -> chi_11 theorem is claimed.
- No legacy physical-role transfer to K_strict_gate is claimed.
- No QW-2191 selector discharge is claimed.
- No ToE closure is claimed.
