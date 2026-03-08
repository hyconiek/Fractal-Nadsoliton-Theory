# N238 Current First Source Topology Nonzero-Flow Subtarget Theorem

Status: `N238_DISCHARGED_CURRENT_FIRST_SOURCE_TOPOLOGY_NONZERO_FLOW_SUBTARGET_THEOREM_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

`N238` packages the current `F130/P218` result into one theorem-level
current-state statement, without upgrading it into:

1. an actual nonzero-flow invariant theorem,
2. a full source-topology nontriviality theorem,
3. a selector-promotion theorem,
4. a `QW-2191` discharge.

## Fixed theorem statement

```text
N238_Current_First_SourceTopology_NonzeroFlow_Subtarget_Theorem

On the current repo state, one explicit future-only nonzero-flow subtarget
packet is exported:

  Xi_src_nonzero_flow_target_v1 : tau_src_candidate_v1 -> source_limit_nonzero_flow_class_v1

where:
  tau_src_candidate_v1
    = (source_limit_tag_v1, phi_barrier_tag_v1, T_flow^(0))

This export is source-side only, observer-free in the witness domain, and
strictly below:
  - actual nonzero-flow discharge,
  - barrier-protected sign discharge,
  - full source-topology nontriviality discharge,
  - selector promotion discharge,
  - quotient-safe QW-2191 resolution,
  - current selector closure,
  - current global QW-2191 discharge.
```

## Why this is the honest theorem

Because on the current repo state:

1. `F127/P215/N235` export only a source-topology candidate packet,
2. `F129/P217/N237` export only a future nontriviality target packet,
3. the next honest unresolved blocker is still the first component:
   actual nonzero-flow of the source packet itself,
4. observer-side asymmetry remains downstream witness only.

Therefore the strongest honest theorem is only:

```text
the current repo exports one explicit future-only nonzero-flow subtarget
for tau_src_candidate_v1
```

and nothing stronger.

## Consequence

After `N238`, the future-route frontier is narrower:

1. the route no longer begins with an unspecified “maybe-nonzero” source class,
2. it begins with one explicit nonzero-flow subtarget packet,
3. but actual nonzero-flow remains open,
4. and all sign, nontriviality, selector-promotion, and `QW-2191` claims remain
   downstream of that open blocker.

## Hard limits

`N238` does not discharge:

1. actual nonzero-flow of `tau_src_candidate_v1`,
2. barrier-protected sign,
3. full source-topology nontriviality,
4. basis-independent selector promotion,
5. quotient-safe `QW-2191` resolution,
6. current selector closure,
7. current global `QW-2191` discharge,
8. ToE closure.
