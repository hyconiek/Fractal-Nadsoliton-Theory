# R25 Direct M2 Shift Equivariance Pairwise Matching Sufficient Route Packet

Status: `R25_EXECUTED_DIRECT_M2_SHIFT_EQUIVARIANCE_PAIRWISE_MATCHING_SUFFICIENT_ROUTE_PACKET_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `R24/P31/N34`, the narrowest route-scoped local blocker on the direct
`m2` branch was:

```text
explicit declared +3 shift-equivariance witness for direct mass-like m2 family
positive support sum
```

`R25` does not pretend to prove that witness.

It attacks only the next honest subobject:

```text
can the direct m2 shift-equivariance target be reduced to a narrower,
explicitly only sufficient, pairwise-matching route on the four direct m2
source-image pairs?
```

This keeps the light issue explicit:

1. the shared kernel/light-facing channel remains already closed by `R14`,
2. `R25` touches only the already exported non-light direct `m2` positive
   support.

## Inputs reused

1. `R24`
   - explicit declared `+3` shift packet on the direct `m2` positive support,
   - explicit source-image coefficient relabeling for the direct `m2` terms.

## Result of `R25`

`R25` materializes one route-scoped sufficient route:

1. `m2_psi1 = m2_psi4`,
2. `m2_psi7 = m2_psi10`,
3. `m2_psi2 = m2_psi5`,
4. `m2_psi8 = m2_psi11`.

If these four pairwise equalities hold, then the direct declared
`+3` shift-equivariance witness holds on the exported positive support sum.

This is only a sufficient route. It is **not** claimed to be necessary or
equivalent.

## Honest frontier after `R25`

On the direct formal coefficient-family route, the host route is reduced to:

1. explicit zero witness for direct quartic-like `g4` family `c1s1`
   shift defect,
2. explicit zero witness for direct quintic-like `g6` family `c1s1`
   shift defect,
3. explicit zero witness for direct yukawa-like `gY` family `c1s1`
   shift defect,
4. four direct `m2` pairwise matching witnesses on one sufficient route,
5. explicit zero witness for the declared `pair1` `c1c1` equation,
6. explicit zero witness for the declared `pair1` `s1s1` equation,
7. `QW-2191` physical canonicalization boundary.

The already closed light-facing kernel channel remains unchanged.

## What `R25` does not claim

`R25` does not claim:

- theorem-level PASS,
- full-closure PASS,
- that any direct `m2` pairwise equality holds,
- that the direct `m2` shift-equivariance holds,
- that the four pairwise conditions are necessary,
- that the sufficient route is equivalent to the direct `m2` shift-equivariance
  target,
- that any other direct family defect vanishes,
- that this direct route is globally equivalent to the main host route,
- that the declared `pair1` `c1c1` equation holds,
- that the declared `pair1` `s1s1` equation holds,
- that `QW-2191` is discharged,
- that selector closure is obtained,
- that ToE is closed.

## Recommended next move

The correct next move is now:

1. rerun the direct formal family route after this route-scoped pairwise
   sufficient packet,
2. accept only:
   - a shorter direct-route missing-object list,
   - or the unchanged negative route.
