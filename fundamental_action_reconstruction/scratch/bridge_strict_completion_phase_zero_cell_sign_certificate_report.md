# Strict completion phase-zero cell-sign certificate

Status: `phase-sign-pattern-derived-from-rational-cell-partition-and-left-anchor`

## Result

This audit derives the integer-node phase sign pattern from the ordered
rational zero-carrier partition and a left anchor sign.  It uses only
rational carrier comparisons and crossing parity, not fresh trigonometric
evaluation of cosine values.

## Sign-rule summary

- Zero carrier order: `['legacy_z0', 'legacy_z1', 'strict_z0', 'legacy_z2']`
- Left anchor sign at d=0: `1`
- Uses trig evaluation: `False`
- All nodes outside zero carriers: `True`
- All node signs match expected: `True`
- All edge flips match crossing parity: `True`
- Maximum crossed zero-carriers per integer edge: `1`
- Derived sign pattern: `[1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1]`
- Derived flip edges: `['1->2', '5->6', '7->8', '9->10']`

## Node sign rows

| d | zero carriers left | parity | derived sign | expected sign | matches |
|---:|---|---|---:|---:|---:|
| 0 | [] | even | 1 | 1 | True |
| 1 | [] | even | 1 | 1 | True |
| 2 | ['legacy_z0'] | odd | -1 | -1 | True |
| 3 | ['legacy_z0'] | odd | -1 | -1 | True |
| 4 | ['legacy_z0'] | odd | -1 | -1 | True |
| 5 | ['legacy_z0'] | odd | -1 | -1 | True |
| 6 | ['legacy_z0', 'legacy_z1'] | even | 1 | 1 | True |
| 7 | ['legacy_z0', 'legacy_z1'] | even | 1 | 1 | True |
| 8 | ['legacy_z0', 'legacy_z1', 'strict_z0'] | odd | -1 | -1 | True |
| 9 | ['legacy_z0', 'legacy_z1', 'strict_z0'] | odd | -1 | -1 | True |
| 10 | ['legacy_z0', 'legacy_z1', 'strict_z0', 'legacy_z2'] | even | 1 | 1 | True |
| 11 | ['legacy_z0', 'legacy_z1', 'strict_z0', 'legacy_z2'] | even | 1 | 1 | True |

## Edge crossing rows

| edge | crossed carriers | count | sign flip | parity match |
|---|---|---:|---:|---:|
| 0->1 | [] | 0 | False | True |
| 1->2 | ['legacy_z0'] | 1 | True | True |
| 2->3 | [] | 0 | False | True |
| 3->4 | [] | 0 | False | True |
| 4->5 | [] | 0 | False | True |
| 5->6 | ['legacy_z1'] | 1 | True | True |
| 6->7 | [] | 0 | False | True |
| 7->8 | ['strict_z0'] | 1 | True | True |
| 8->9 | [] | 0 | False | True |
| 9->10 | ['legacy_z2'] | 1 | True | True |
| 10->11 | [] | 0 | False | True |

## Proof certificate

- `anchor_step`: The left anchor sign at d=0 is fixed as +1 from the already-certified phase pattern context.
- `counting_step`: For each integer node, count ordered zero-carriers strictly to the left using rational interval comparisons only.
- `node_step`: Node-clearance guarantees no integer node lies inside a zero-carrier interval, so each node has a well-defined cell sign.
- `edge_step`: Adjacent node signs differ exactly when an odd number of zero-carriers is crossed by the edge.
- `nonduplication`: This is a combinatorial cell-sign theorem, not another zero-location, node-clearance, cell-partition, damping, cocycle, or chain-integrity audit.
- `theoretical_limit`: The certificate derives the finite sign pattern from the selected zero-carrier partition; it does not derive omega/phi or transport from strict nadsoliton dynamics.

## Hard limits

- K_strict_gate remains the current live/full operational kernel.
- No unqualified identity K_legacy_ont == K_strict_gate is claimed.
- No proof derives strict omega/phi or phase transport from strict nadsoliton dynamics.
- No beta_tors -> chi_11 theorem is claimed.
- No legacy physical-role transfer to K_strict_gate is claimed.
- No QW-2191 selector discharge is claimed.
- No ToE closure is claimed.
