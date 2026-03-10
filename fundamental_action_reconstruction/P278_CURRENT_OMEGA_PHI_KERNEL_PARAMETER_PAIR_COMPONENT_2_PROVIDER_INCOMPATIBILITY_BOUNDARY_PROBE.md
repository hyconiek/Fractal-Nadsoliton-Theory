# P278 Current Omega-Phi Kernel Parameter Pair Component 2 Provider Incompatibility Boundary Probe

Status: `P278_EXECUTED_CURRENT_OMEGA_PHI_KERNEL_PARAMETER_PAIR_COMPONENT_2_PROVIDER_INCOMPATIBILITY_BOUNDARY_PROBE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

Probe whether the current repo may already treat `(omega,phi)` as the
pair-indexed population anchor for component 2.

## Probe table

| Check | Verdict | Meaning |
|---|---|---|
| `(omega,phi)` exist as strict-kernel parameters | YES | the pair is real |
| `(omega,phi)` support local source-side derivative calculus | YES | `F163/N285` already use them there |
| `(omega,phi)` are typed as pair-indexed carrier over `[pair1,pair2]` | NO | no such typing is exported |
| `(omega,phi)` export a typed map to `theta_1/theta_2` | NO | no such map is exported |
| `(omega,phi)` export a basis-pair population rule | NO | no such rule is exported |
| `(omega,phi)` export actual component-2 support | NO | the proposal remains below entry |

## Probe result

`P278` returns:

```text
omega_phi_pair_exists: yes
omega_phi_component_2_anchor_support: no
```

## Consequence

The strongest honest current repo reading is:

1. `(omega,phi)` may be cited as a real source-side parameter pair,
2. but not yet as the pair-indexed population anchor required by component 2.
