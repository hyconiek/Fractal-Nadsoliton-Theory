# N247 Current First Actual Source Topology Barrier Sign Component Witness Theorem

Status: `N247_DISCHARGED_CURRENT_FIRST_ACTUAL_SOURCE_TOPOLOGY_BARRIER_SIGN_COMPONENT_WITNESS_THEOREM_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

`N247` packages the current `F139/P227` result into one theorem-level
current-state statement, without upgrading it into:

1. a full barrier-protected sign theorem,
2. a full source-topology nontriviality theorem,
3. a selector-promotion theorem,
4. a quotient-safe `QW-2191` discharge.

## Fixed theorem statement

```text
N247_Current_First_Actual_SourceTopology_BarrierSign_Component_Witness_Theorem

On the current repo state, one explicit scalar barrier margin and one actual
source-side scalar barrier-sign component witness are exported:

  delta_src_barrier_sign_margin_v1 := (pi/2) - |phi| = 1.4082963267948965 > 0
  psi_src_barrier_sign_component_witness_v1 := sign(cos(phi)) = +1

for the candidate packet

  tau_src_candidate_v1
    = (source_limit_tag_v1, phi_barrier_tag_v1, T_flow^(0))

This witness remains:
  - observer-free in the witness domain,
  - below full barrier-protected sign discharge,
  - below full source-topology nontriviality discharge,
  - below basis-independent selector promotion,
  - below quotient-safe QW-2191 resolution,
  - below current selector closure,
  - below current global QW-2191 discharge.
```

## Why this is the honest theorem

Because on the current repo state:

1. `F131/P219/N239` export only a future barrier-sign subtarget,
2. `F138/P226/N246` add one actual scalar nonzero-flow component witness,
3. `F139/P227` add one actual scalar sign component witness together with one
   explicit positive barrier margin on the current declared core branch,
4. but full barrier-protected sign discharge, full source-topology
   nontriviality, basis independence, selector promotion, and `QW-2191` safety
   remain open,
5. observer-side asymmetry remains downstream witness only.

Therefore the strongest honest theorem is only:

```text
the current repo exports one actual scalar barrier-sign component witness
for tau_src_candidate_v1
```

and nothing stronger.

## Consequence

After `N247`, the future-route frontier is narrower:

1. the route no longer contains only an abstract sign subtarget below
   `tau_src_candidate_v1`,
2. it now contains one actual scalar sign component witness with an explicit
   positive barrier margin on the current core branch,
3. but full barrier-protected sign discharge, full source-topology
   nontriviality, selector promotion, and `QW-2191` resolution remain
   downstream of that still-limited witness.

## Hard limits

`N247` does not discharge:

1. full barrier-protected sign,
2. full source-topology nontriviality,
3. basis-independent selector promotion,
4. quotient-safe `QW-2191` resolution,
5. current selector closure,
6. current global `QW-2191` discharge,
7. ToE closure.
