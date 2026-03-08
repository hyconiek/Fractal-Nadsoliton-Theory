# RELEASE 5.4: First Actual Source-Topology Component Witness

**Version:** 5.4.0  
**Date:** 2026-03-08  
**Branch:** `main`

## Executive Summary

Release 5.4 does not close the Theory of Everything and does not discharge
`QW-2191`.

What it adds is narrower and scientifically cleaner: after Release 5.3 opened
the `T14` Source Topology Selector route only as a future theorem-spec, Release
5.4 exports the first actual source-side component witness on that route.

The key new result is:

```text
xi_src_nonzero_flow_component_witness_v1 := |cos(phi)| = 0.9868259031903286 > 0
```

This is still strictly below:

1. barrier-protected sign discharge,
2. full source-topology nontriviality,
3. basis-independent selector promotion,
4. quotient-safe `QW-2191` resolution,
5. global selector closure,
6. ToE closure.

So Release 5.4 is not a release of final closure. It is the release where the
`T14` lane ceases to be purely abstract and acquires its first actual
source-side witness.

## 1. What Changes Relative to Release 5.3

Release 5.3 established:

1. the observer-side global-promotion obstruction in
   `N234`,
2. the `T14` Source Topology Selector route as a future theorem-spec,
3. a decomposition of the missing source-topology ingredients into:
   nonzero flow, barrier-protected sign, observer-free scope, basis
   independence, and quotient-safe promotion.

Release 5.4 adds the first actual component-level step below that decomposition.

The route now has:

1. one explicit source-topology candidate packet,
2. one explicit selector-promotion target,
3. one explicit nontriviality-component package,
4. one actual source-side scalar witness for nonzero flow.

This is a real methodological shift:

- Release 5.3 said what must eventually be proved,
- Release 5.4 proves one first component-level fact that sits strictly inside
  that future route.

## 2. Core Result: The First Actual Source-Topology Witness

### 2.1 Candidate packet

The current source-topology candidate packet remains:

```text
tau_src_candidate_v1 :=
(
  source_limit_tag_v1,
  phi_barrier_tag_v1,
  T_flow^(0)
)
```

with

```text
K_strict(0) = cos(phi) = 0.9868259031903286
T_flow^(0) = cos(phi) * e_topo
```

### 2.2 First actual scalar component witness

Release 5.4 exports:

```text
xi_src_nonzero_flow_component_witness_v1 := |cos(phi)|
```

Numerically:

```text
xi_src_nonzero_flow_component_witness_v1 = 0.9868259031903286 > 0
```

This is the first actual source-side scalar witness on the `T14` route.

It means the route no longer begins only with an abstract nonzero-flow class.
It now contains one explicit current witness that the source-limit packet
contains a non-vanishing source-side component.

## 3. Why This Matters

### 3.1 It is upstream of the observer

The new witness is extracted before all observer-side pushforward and readout
layers.

So the preferred causal order remains unchanged:

```text
nadsoliton -> light -> matter -> emergent observer
```

The observer is still downstream only.

### 3.2 It does not rely on observer promotion

The witness is not inferred from observer stability, readout, closure-candidate
chains, or downstream coarse graining.

This is important because `N234` already blocks the false move:

```text
stable downstream observer asymmetry
!= theorem-level global selector closure
```

Release 5.4 therefore strengthens the project exactly where it must be
strengthened: on the source side, not on the observer side.

### 3.3 It is still kernel-split-safe

The witness uses the strict source-limit control datum only as currently
exported. It does not:

1. identify `K_strict_gate` with `K_legacy_ont`,
2. transfer legacy physical-role claims,
3. claim a legacy-to-strict bridge.

So it remains compatible with `K1`, `K2`, `F2`, `F3`, and `S2`.

## 4. What Release 5.4 Proves

Release 5.4 proves, on the current repo state, only the following stronger but
still limited statement:

1. the `T14` route now contains one actual source-side scalar nonzero-flow
   component witness,
2. that witness is observer-free in its witness domain,
3. the witness remains below all theorem-level promotions that would be needed
   for selector closure or `QW-2191` discharge.

The theorem-level packaging of this result is:

- `F138`
- `P226`
- `N246`

## 5. What Release 5.4 Does Not Prove

Release 5.4 still does not prove:

1. actual barrier-protected sign discharge,
2. full source-topology nontriviality,
3. basis-independent selector promotion,
4. quotient-safe `QW-2191` resolution,
5. current global `QW-2191` discharge,
6. strict-core selector closure,
7. final actual emergent-observer closure,
8. final legacy-to-strict bridge,
9. final ToE closure.

It also does not promote the observer into the source of asymmetry.

Current reading remains:

- source-side topology is where the honest reopening route lives,
- observer-side asymmetry remains downstream structural evidence only.

## 6. Exact Next Step

The exact next honest move after Release 5.4 is:

1. attempt an actual barrier-protected sign witness on the `T14` route,
2. only after that attempt a basis-independent `Pi_sel` promotion witness,
3. only after both attempt any quotient-safe `QW-2191` resolution.

The project should not jump directly from the scalar witness to global selector
closure.

## 7. Main Artifacts

- `fundamental_action_reconstruction/README.md`
- `N234_CURRENT_GLOBAL_SELECTOR_CLOSURE_AND_QW2191_DISCHARGE_PROMOTION_OBSTRUCTION_THEOREM.md`
- `T14_SOURCE_TOPOLOGY_SELECTOR_THEOREM_SPEC.md`
- `F127_FIRST_SOURCE_TOPOLOGY_INVARIANT_CANDIDATE_PACKET.md`
- `N235_CURRENT_FIRST_SOURCE_TOPOLOGY_INVARIANT_CANDIDATE_PACKET_THEOREM.md`
- `F128_FIRST_SOURCE_TOPOLOGY_SELECTOR_PROMOTION_TARGET_PACKET.md`
- `N236_CURRENT_FIRST_SOURCE_TOPOLOGY_SELECTOR_PROMOTION_TARGET_THEOREM.md`
- `F133_FIRST_SOURCE_TOPOLOGY_NONTRIVIALITY_COMPONENTS_PACKAGE_PACKET.md`
- `N241_CURRENT_FIRST_SOURCE_TOPOLOGY_NONTRIVIALITY_COMPONENTS_PACKAGE_THEOREM.md`
- `F137_FIRST_SOURCE_TOPOLOGY_QUOTIENT_SAFE_QW2191_RESOLUTION_TARGET_PACKET.md`
- `N245_CURRENT_FIRST_SOURCE_TOPOLOGY_QUOTIENT_SAFE_QW2191_RESOLUTION_TARGET_THEOREM.md`
- `F138_FIRST_ACTUAL_SOURCE_TOPOLOGY_NONZERO_FLOW_COMPONENT_WITNESS_PACKET.md`
- `P226_CURRENT_ACTUAL_SOURCE_TOPOLOGY_NONZERO_FLOW_COMPONENT_WITNESS_PROBE.md`
- `N246_CURRENT_FIRST_ACTUAL_SOURCE_TOPOLOGY_NONZERO_FLOW_COMPONENT_WITNESS_THEOREM.md`
