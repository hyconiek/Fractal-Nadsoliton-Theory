# N251 Current First Actual Source Topology Nonzero-Flow Witness Theorem

Status: `N251_DISCHARGED_CURRENT_FIRST_ACTUAL_SOURCE_TOPOLOGY_NONZERO_FLOW_WITNESS_THEOREM_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

`N251` packages the current `F143/P231` result into one theorem-level
current-state statement, without upgrading it into:

1. a full source-topology nontriviality theorem,
2. a selector-promotion theorem,
3. a quotient-safe `QW-2191` discharge.

## Fixed theorem statement

```text
N251_Current_First_Actual_SourceTopology_NonzeroFlow_Witness_Theorem

On the current repo state, one actual source-side nonzero-flow witness is
exported:

  Xi_src_nonzero_flow_actual_witness_v1 :
    tau_src_candidate_v1 -> source_limit_nonzero_flow_class_v1

supported by the current-repo-state packet

  W_src_nonzero_flow_support_packet_v1
    = (
        source_limit_tag_v1,
        T_flow^(0),
        xi_src_nonzero_flow_component_witness_v1
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

1. `F130/P218/N238` export only a future nonzero-flow subtarget,
2. `F138/P226/N246` add one actual scalar nonzero-flow component witness,
3. `F143/P231` add only the packet-level lift from that already exported
   scalar witness into the already declared class
   `source_limit_nonzero_flow_class_v1`,
4. but full source-topology nontriviality, basis independence, selector
   promotion, and `QW-2191` safety remain open,
5. observer-side asymmetry remains downstream witness only.

Therefore the strongest honest theorem is:

```text
the current repo exports one actual source-side nonzero-flow witness
for tau_src_candidate_v1
```

and nothing stronger.

## Consequence

After `N251`, the future-route frontier is narrower:

1. the nonzero-flow layer no longer remains only a future subtarget or a
   merely scalar component witness,
2. it now contains one actual witness in
   `source_limit_nonzero_flow_class_v1`,
3. but full source-topology nontriviality, basis-independent selector
   promotion, and `QW-2191` resolution remain downstream of that still-limited
   witness.

## Hard limits

`N251` does not discharge:

1. full source-topology nontriviality,
2. basis-independent selector promotion,
3. quotient-safe `QW-2191` resolution,
4. current selector closure,
5. current global `QW-2191` discharge,
6. ToE closure.
