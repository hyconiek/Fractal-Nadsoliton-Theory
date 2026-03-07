# R26 Direct M2 Psi1 Psi4 Role Matching Packet

Status: `R26_EXECUTED_DIRECT_M2_PSI1_PSI4_ROLE_MATCHING_PACKET_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `R25/P32/N35`, the narrowest route-scoped local blocker on the direct
`m2` sufficient lane included the single pairwise target:

```text
explicit pairwise matching witness for m2_psi1 = m2_psi4
```

`R26` does not pretend to prove that equality.

It attacks only the next honest subobject:

```text
can the single pairwise target m2_psi1 = m2_psi4 be reduced to one declared
action/eom role-matching packet under the already exported plus3 shift?
```

This keeps the light issue explicit:

1. the shared kernel/light-facing channel remains already closed by `R14`,
2. `R26` touches only one non-light direct `m2` pair on the already exported
   sufficient route.

## Inputs reused

1. `R24`
   - explicit declared `psi1 -> psi4`,
   - explicit declared `m2_psi1 -> m2_psi4`.
2. `QW-2164`
   - canonical continuum potential with explicit quadratic `m2` terms.
3. `QW-2165`
   - exhaustive canonical Euler-Lagrange system with explicit local `m2`
     terms in `eom_psi1` and `eom_psi4`.

## Result of `R26`

`R26` materializes one route-scoped declared role-matching packet:

1. source action term:
   `m2_psi1*psi1**2/2`,
2. target action term:
   `m2_psi4*psi4**2/2`,
3. declared shifted source action term:
   `m2_psi4*psi4**2/2`,
4. source local eom term:
   `m2_psi1*psi1(x)`,
5. target local eom term:
   `m2_psi4*psi4(x)`,
6. declared shifted source local eom term:
   `m2_psi4*psi4(x)`.

So, on this one pair only, the repo now exports exact slot-role matching in the
canonical action and in the local eom under the declared `+3` shift.

This is still **not** a coefficient-equality witness.

## Honest frontier after `R26`

On the direct formal coefficient-family route, the host route is reduced to:

1. explicit zero witness for direct quartic-like `g4` family `c1s1`
   shift defect,
2. explicit zero witness for direct quintic-like `g6` family `c1s1`
   shift defect,
3. explicit zero witness for direct yukawa-like `gY` family `c1s1`
   shift defect,
4. explicit coefficient-identification witness for the declared role-matched
   direct `m2` pair `m2_psi1 = m2_psi4`,
5. explicit pairwise matching witness for `m2_psi7 = m2_psi10`,
6. explicit pairwise matching witness for `m2_psi2 = m2_psi5`,
7. explicit pairwise matching witness for `m2_psi8 = m2_psi11`,
8. explicit zero witness for the declared `pair1` `c1c1` equation,
9. explicit zero witness for the declared `pair1` `s1s1` equation,
10. `QW-2191` physical canonicalization boundary.

The already closed light-facing kernel channel remains unchanged.

## What `R26` does not claim

`R26` does not claim:

- theorem-level PASS,
- full-closure PASS,
- that `m2_psi1 = m2_psi4`,
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

1. rerun the direct formal family route after this one-pair role-matching
   packet,
2. accept only:
   - a shorter direct-route missing-object list,
   - or the unchanged negative route.
