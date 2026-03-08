# N249 Current First Actual Source Topology Barrier-Protected Sign Witness Theorem

Status: `N249_DISCHARGED_CURRENT_FIRST_ACTUAL_SOURCE_TOPOLOGY_BARRIER_PROTECTED_SIGN_WITNESS_THEOREM_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

`N249` packages the current `F141/P229` result into one theorem-level
current-state statement, without upgrading it into:

1. a full source-topology nontriviality theorem,
2. a selector-promotion theorem,
3. a quotient-safe `QW-2191` discharge.

## Fixed theorem statement

```text
N249_Current_First_Actual_SourceTopology_BarrierProtectedSign_Witness_Theorem

On the current repo state, one actual source-side barrier-protected sign
witness is exported:

  Psi_src_barrier_sign_actual_witness_v1 :
    tau_src_candidate_v1 -> barrier_protected_sign_class_v1

supported by the current-repo-state packet

  W_src_barrier_sign_support_packet_v1
    = (
        phi_barrier_tag_v1,
        delta_src_barrier_sign_margin_v1,
        epsilon_src_local_barrier_radius_v1,
        psi_src_barrier_sign_component_witness_v1,
        chi_src_local_barrier_sign_stability_witness_v1
      )

This witness remains:
  - observer-free in the witness domain,
  - below full source-topology nontriviality discharge,
  - below basis-independent selector promotion,
  - below quotient-safe QW-2191 resolution,
  - below current selector closure,
  - below current global QW-2191 discharge.
```

## Why this is the honest theorem

Because on the current repo state:

1. `F131/P219/N239` export only a future sign subtarget,
2. `F139/P227/N247` add actual scalar barrier-sign support,
3. `F140/P228/N248` add actual positive-radius local sign stability,
4. `F141/P229` add only the packet-level lift from those already exported
   support objects into the already declared class
   `barrier_protected_sign_class_v1`,
5. but full source-topology nontriviality, basis independence, selector
   promotion, and `QW-2191` safety remain open,
6. observer-side asymmetry remains downstream witness only.

Therefore the strongest honest theorem is:

```text
the current repo exports one actual source-side barrier-protected sign witness
for tau_src_candidate_v1
```

and nothing stronger.

## Consequence

After `N249`, the future-route frontier is narrower:

1. the sign layer no longer remains only a future subtarget or a merely local
   branch witness,
2. it now contains one actual witness in
   `barrier_protected_sign_class_v1`,
3. but full source-topology nontriviality, basis-independent selector
   promotion, and `QW-2191` resolution remain downstream of that still-limited
   witness.

## Hard limits

`N249` does not discharge:

1. full source-topology nontriviality,
2. basis-independent selector promotion,
3. quotient-safe `QW-2191` resolution,
4. current selector closure,
5. current global `QW-2191` discharge,
6. ToE closure.
