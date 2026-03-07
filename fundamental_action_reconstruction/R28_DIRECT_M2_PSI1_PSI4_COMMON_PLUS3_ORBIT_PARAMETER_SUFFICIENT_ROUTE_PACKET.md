# R28 Direct M2 Psi1 Psi4 Common Plus3 Orbit Parameter Sufficient Route Packet

Status: `R28_EXECUTED_DIRECT_M2_PSI1_PSI4_COMMON_PLUS3_ORBIT_PARAMETER_SUFFICIENT_ROUTE_PACKET_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `R27/P34/N37`, the narrowest route-scoped local blocker on the single
direct `m2` pair was:

```text
explicit common-parameter-source or symbol-identification witness for the
declared formal m2 slots m2_psi1 and m2_psi4
```

`R28` does not pretend to prove that witness.

It attacks only the next honest subobject:

```text
can that one common-source/symbol-identification gap be reduced to one even
narrower sufficient route through a single common plus3 carrier-segment
parameter for m2_psi1 / m2_psi4?
```

This keeps the light issue explicit:

1. the shared kernel/light-facing channel remains already closed by `R14`,
2. `R28` touches only one non-light direct `m2` pair on the already exported
   sufficient route.

## Inputs reused

1. `R24`
   - declared `psi1 -> psi4`,
   - declared `m2_psi1 -> m2_psi4`.
2. `R27`
   - declared formal slot separation for `m2_psi1 / m2_psi4`.

## Result of `R28`

`R28` materializes one route-scoped sufficient route:

1. `m2_psi1 = mu_m2_plus3_segment_psi1_psi4`,
2. `m2_psi4 = mu_m2_plus3_segment_psi1_psi4`.

If these two assignments hold, then the pairwise equality
`m2_psi1 = m2_psi4` follows on this route.

This is only a sufficient route. It is **not** claimed to be necessary or
equivalent.

## Honest frontier after `R28`

On the direct formal coefficient-family route, the host route is reduced to:

1. explicit zero witness for direct quartic-like `g4` family `c1s1`
   shift defect,
2. explicit zero witness for direct quintic-like `g6` family `c1s1`
   shift defect,
3. explicit zero witness for direct yukawa-like `gY` family `c1s1`
   shift defect,
4. explicit assignment witness of `m2_psi1` and `m2_psi4` to one common
   plus3 carrier-segment parameter,
5. explicit pairwise matching witness for `m2_psi7 = m2_psi10`,
6. explicit pairwise matching witness for `m2_psi2 = m2_psi5`,
7. explicit pairwise matching witness for `m2_psi8 = m2_psi11`,
8. explicit zero witness for the declared `pair1` `c1c1` equation,
9. explicit zero witness for the declared `pair1` `s1s1` equation,
10. `QW-2191` physical canonicalization boundary.

The already closed light-facing kernel channel remains unchanged.

## What `R28` does not claim

`R28` does not claim:

- theorem-level PASS,
- full-closure PASS,
- that `m2_psi1 = m2_psi4`,
- that a common plus3 carrier-segment parameter actually exists,
- that the sufficient route is necessary or equivalent,
- that any other direct `m2` pairwise equality holds,
- that the direct `m2` shift-equivariance holds,
- that any direct `g4/g6/gY` family defect vanishes,
- that this direct route is globally equivalent to the main host route,
- that the declared `pair1` `c1c1` equation holds,
- that the declared `pair1` `s1s1` equation holds,
- that `QW-2191` is discharged,
- that selector closure is obtained,
- that ToE is closed.

## Recommended next move

The correct next move is now:

1. rerun the direct formal family route after this one-pair common-parameter
   sufficient packet,
2. accept only:
   - a shorter direct-route missing-object list,
   - or the unchanged negative route.
