# R31 Direct M2 Psi1 Source Action Common Monomial Support Packet

Status: `R31_EXECUTED_DIRECT_M2_PSI1_SOURCE_ACTION_COMMON_MONOMIAL_SUPPORT_PACKET_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `R30/P37/N40`, the narrowest route-scoped blocker on the attacked source
action side was:

```text
explicit source action-role assignment witness for m2_psi1 to
mu_m2_plus3_segment_psi1_psi4
```

`R31` does not pretend to prove that source action-role assignment witness.

It attacks only the next honest subobject:

```text
can that one source action-role assignment witness be reduced to one
coefficient-identification witness on the already fixed common source-action
monomial support psi1**2/2?
```

This keeps the light issue explicit:

1. the shared kernel/light-facing channel remains already closed by `R14`,
2. `R31` touches only one non-light direct `m2` source action-role on the
   already exported sufficient route.

## Inputs reused

1. `R30`
   - exact source action-role assignment witness under attack,
   - declared lifted source action term.
2. `R26`
   - exact source action term `m2_psi1*psi1**2/2`.

## Result of `R31`

`R31` materializes one exact route-scoped common-monomial-support packet:

1. source action term:
   `m2_psi1*psi1**2/2`,
2. declared lifted action term:
   `mu_m2_plus3_segment_psi1_psi4*psi1**2/2`,
3. fixed common monomial support:
   `psi1**2/2`,
4. remaining missing witness:
   `explicit_source_action_monomial_coefficient_identification_witness_for_m2_psi1_and_mu_m2_plus3_segment_psi1_psi4_on_common_psi1_squared_over_2_support`.

So, on this one source action-role only, the repo now exports not only the
role-specific assignment gap but also the exact current reason why it still
fails: the coefficient labels on the already fixed monomial support are still
not identified.

This is **not** a global cancellation argument and **not** a nonzero-factor
claim.

## Honest frontier after `R31`

On the direct formal coefficient-family route, the host route is reduced to:

1. explicit zero witness for direct quartic-like `g4` family `c1s1`
   shift defect,
2. explicit zero witness for direct quintic-like `g6` family `c1s1`
   shift defect,
3. explicit zero witness for direct yukawa-like `gY` family `c1s1`
   shift defect,
4. explicit source action monomial coefficient-identification witness for
   `m2_psi1` and `mu_m2_plus3_segment_psi1_psi4` on common
   `psi1**2/2` support,
5. explicit source eom-role assignment witness for `m2_psi1` to
   `mu_m2_plus3_segment_psi1_psi4`,
6. explicit assignment witness of `m2_psi4` to
   `mu_m2_plus3_segment_psi1_psi4`,
7. explicit pairwise matching witness for `m2_psi7 = m2_psi10`,
8. explicit pairwise matching witness for `m2_psi2 = m2_psi5`,
9. explicit pairwise matching witness for `m2_psi8 = m2_psi11`,
10. explicit zero witness for the declared `pair1` `c1c1` equation,
11. explicit zero witness for the declared `pair1` `s1s1` equation,
12. `QW-2191` physical canonicalization boundary.

The already closed light-facing kernel channel remains unchanged.

## What `R31` does not claim

`R31` does not claim:

- theorem-level PASS,
- full-closure PASS,
- that `m2_psi1 = m2_psi4`,
- that any common plus3 carrier-segment parameter actually exists,
- that the source action-side coefficient-identification witness is present,
- that any global cancellation or nonzero-factor argument holds,
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

1. rerun the direct formal family route after this exact common-support packet,
2. accept only:
   - a shorter direct-route missing-object list,
   - or the unchanged negative route.
