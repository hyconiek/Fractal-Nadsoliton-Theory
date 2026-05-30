# Strict completion phase-zero cell-partition certificate

Status: `phase-zero-carriers-form-positive-rational-cell-partition-on-audited-domain`

## Result

This audit proves that the in-domain phase-zero carriers are ordered,
pairwise separated, and strictly inside open integer edges.  The resulting
cell partition recovers the established phase flip edges and sign pattern,
without claiming a strict dynamical derivation of the phase parameters.

## Cell-partition summary

- Domain zero carrier order: `['legacy_z0', 'legacy_z1', 'strict_z0', 'legacy_z2']`
- All carriers inside open edges: `True`
- Adjacent carriers ordered/disjoint: `True`
- All cells positive length: `True`
- Minimum boundary clearance: `1/3` = `3.333333333333e-01`
- Minimum adjacent zero separation: `441202/251877` = `1.751656562529e+00`
- Minimum cell length: `4/3` = `1.333333333333e+00`
- Derived flip edges: `['1->2', '5->6', '7->8', '9->10']`
- Derived sign pattern: `[1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1]`

## Edge zero inventory

| edge | zero carriers | count | odd parity flip |
|---|---|---:|---:|
| 0->1 | [] | 0 | False |
| 1->2 | ['legacy_z0'] | 1 | True |
| 2->3 | [] | 0 | False |
| 3->4 | [] | 0 | False |
| 4->5 | [] | 0 | False |
| 5->6 | ['legacy_z1'] | 1 | True |
| 6->7 | [] | 0 | False |
| 7->8 | ['strict_z0'] | 1 | True |
| 8->9 | [] | 0 | False |
| 9->10 | ['legacy_z2'] | 1 | True |
| 10->11 | [] | 0 | False |

## Adjacent zero separation

| left | right | separation | ordered/disjoint |
|---|---|---:|---:|
| legacy_z0 | legacy_z1 | 4/1 | True |
| legacy_z1 | strict_z0 | 265586/118137 | True |
| strict_z0 | legacy_z2 | 441202/251877 | True |

## Proof certificate

- `carrier_step`: Legacy exact zeros and strict interval zeros are represented as rational zero-carriers on the audited line.
- `edge_step`: Every in-domain zero-carrier is strictly inside a single open integer edge with positive rational boundary clearance.
- `separation_step`: Adjacent in-domain zero-carriers have positive rational separation, so no two phase-zero events collide or overlap.
- `cell_step`: The ordered carriers cut [0,11] into positive-length rational cells; signs can change only across odd-parity carrier edges.
- `parity_step`: The induced odd-parity edge inventory recovers the established flip edges and sign pattern.
- `nonduplication`: This is a cell-partition/order/separation certificate, not another node-clearance, margin, damping, cocycle, or chain-integrity audit.
- `theoretical_limit`: The certificate proves finite rational phase-zero carrier separation; it does not derive omega/phi or transport from strict nadsoliton dynamics.

## Hard limits

- K_strict_gate remains the current live/full operational kernel.
- No unqualified identity K_legacy_ont == K_strict_gate is claimed.
- No proof derives strict omega/phi or phase transport from strict nadsoliton dynamics.
- No beta_tors -> chi_11 theorem is claimed.
- No legacy physical-role transfer to K_strict_gate is claimed.
- No QW-2191 selector discharge is claimed.
- No ToE closure is claimed.
