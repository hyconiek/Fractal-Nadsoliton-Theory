# R10 Legacy Control To Current-Pair Chart Reduction Packet For K_obs

Status: `R10_EXECUTED_LEGACY_CONTROL_TO_CURRENT_PAIR_CHART_REDUCTION_PACKET_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `P13/N16`, the next missing factorization subobject was:

```text
selector_sector_reduction_of_the_legacy_control_side_onto_pair1_or_an_equivalent_actual_target
```

`R10` tests the narrowest honest constructive question:

```text
does the current repo already contain a typed reduction from the legacy
control-side carrier M_control into the chosen explicit current-pair chart
V_1 = span{c1,s1}, even if that chart still lacks strict selector privilege?
```

`R10` does not claim strict selector-target justification.

## Inputs reused

1. `R9`
   - typed host-to-control pushforward
     `P_control : Psi_host_12 -> M_control`.
2. `C15`
   - declared control basis
     `B_control = (c1,s1,c2,s2)`.
3. `C29`
   - explicit local reduced-projector formula on each pair plane.
4. `H8`
   - explicit `H3` chain starts on
     `V_1 = span{c1,s1}`.
5. `H33/H34`
   - `pair1` remains only a local chart and not a strict selector target.

## Result of `R10`

`R10` establishes:

1. the control-side carrier basis is already explicit:
   `(c1,s1,c2,s2)`,
2. the explicit current-pair chain domain is already explicit:
   `V_1 = span{c1,s1}`,
3. the typed chart-reduction map can already be packetized as:

```text
Pi_pair1 : M_control -> V_1
Pi_pair1 = [[1,0,0,0],[0,1,0,0]]
```

with basis order
`(c1,s1,c2,s2) -> (c1,s1)`.

So the reduction blocker is now resolved at **chosen current-pair chart**
level:

```text
legacy control side -> chosen explicit current-pair chart = present
```

## Honest frontier after `R10`

`R10` does **not** establish:

- that `pair1` is a privileged selector target,
- a basis-covariant or target-independent selector reduction,
- an intertwiner/equality witness to the computed `P10` current-pair block,
- factorization of existing kernel feedback into `K_obs`.

So the route frontier becomes:

- `R10_B1 := typed reduction into the chosen explicit current-pair chart is now present, but no intertwiner/equality witness identifies the resulting chart-reduced legacy side with the computed current-pair H3 block`

## What `R10` does not claim

`R10` does not claim:

- theorem-level PASS,
- full-closure PASS,
- that `pair1` is globally or physically canonical,
- that strict selector-sector reduction is discharged,
- that `QW-2191` is discharged,
- that ToE is closed.

## Recommended next move

The correct next move is now:

1. rerun the exact same factorization route after this chart-reduction packet,
2. accept only:
   - a single residual missing-object witness,
   - or the unchanged negative route.
