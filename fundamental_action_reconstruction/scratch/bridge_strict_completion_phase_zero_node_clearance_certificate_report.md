# Strict completion phase-zero node-clearance certificate

Status: `all-integer-phase-nodes-have-positive-rational-clearance-from-legacy-and-strict-zeros`

## Result

This audit gives rational lower bounds proving that no audited integer
node `d=0..11` is accidentally sitting on a legacy or strict phase zero.
It complements edge-parity and parameter-margin certificates without
claiming a strict dynamical derivation of the phase parameters.

## Clearance summary

- All legacy zeros off integer nodes: `True`
- All strict zero intervals off integer nodes: `True`
- All audited integer nodes certified nonzero: `True`
- Minimum legacy clearance: `1/3` = `3.333333333333e-01`
- Minimum strict clearance lower bound: `35122/83959` = `4.183232291952e-01`
- Minimum combined node clearance lower bound: `1/3` = `3.333333333333e-01`
- Preserved flip edges: `['1->2', '5->6', '7->8', '9->10']`
- Preserved sign pattern: `[1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1]`

## Strict zero clearance rows

| k | location | lower | upper | min node clearance lower bound | nearest nodes |
|---:|---|---:|---:|---:|---|
| -1 | left-of-domain | -783450/83959 | -367450/39379 | 367450/39379 | [0] |
| 0 | 7->8 | 298550/39379 | 636550/83959 | 35122/83959 | [8] |
| 1 | right-of-domain | 964550/39379 | 2056550/83959 | 531381/39379 | [11] |

## Integer node rows

| d | legacy clearance | strict clearance lower bound | combined lower bound | certified nonzero |
|---:|---:|---:|---:|---:|
| 0 | 4/3 | 298550/39379 | 4/3 | True |
| 1 | 1/3 | 259171/39379 | 1/3 | True |
| 2 | 2/3 | 219792/39379 | 2/3 | True |
| 3 | 5/3 | 180413/39379 | 5/3 | True |
| 4 | 4/3 | 141034/39379 | 4/3 | True |
| 5 | 1/3 | 101655/39379 | 1/3 | True |
| 6 | 2/3 | 62276/39379 | 2/3 | True |
| 7 | 5/3 | 22897/39379 | 22897/39379 | True |
| 8 | 4/3 | 35122/83959 | 35122/83959 | True |
| 9 | 1/3 | 119081/83959 | 1/3 | True |
| 10 | 2/3 | 203040/83959 | 2/3 | True |
| 11 | 5/3 | 286999/83959 | 5/3 | True |

## Proof certificate

- `legacy_clearance_step`: Legacy zeros are exact rationals 4/3, 16/3, 28/3, hence their minimum integer-node clearance is exactly positive.
- `strict_clearance_step`: Each relevant strict zero is enclosed using 333/106 < pi < 355/113 and omega=743/4000, phi=13/80; integer-node distance to the whole interval gives a rational lower bound.
- `node_step`: For each d=0..11, the minimum of all legacy and strict zero clearances is positive, so no audited integer node is a phase zero.
- `sign_step`: With damping already positive, sign changes in the completion factor remain controlled by phase-zero parity rather than damping or node degeneracy.
- `nonduplication`: This is a node-clearance/no-degeneracy certificate, not another parameter-robustness, damping-monotonicity, transport-cocycle, or chain-integrity audit.
- `theoretical_limit`: The certificate proves finite rational separation of selected phase zeros from integer nodes; it does not derive omega/phi or transport from strict nadsoliton dynamics.

## Hard limits

- K_strict_gate remains the current live/full operational kernel.
- No unqualified identity K_legacy_ont == K_strict_gate is claimed.
- No proof derives strict omega/phi or phase transport from strict nadsoliton dynamics.
- No beta_tors -> chi_11 theorem is claimed.
- No legacy physical-role transfer to K_strict_gate is claimed.
- No QW-2191 selector discharge is claimed.
- No ToE closure is claimed.
