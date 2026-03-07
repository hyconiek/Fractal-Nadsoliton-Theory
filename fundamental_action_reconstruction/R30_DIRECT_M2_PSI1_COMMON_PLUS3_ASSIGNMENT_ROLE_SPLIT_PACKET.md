# R30 Direct M2 Psi1 Common Plus3 Assignment Role Split Packet

Status: `R30_EXECUTED_DIRECT_M2_PSI1_COMMON_PLUS3_ASSIGNMENT_ROLE_SPLIT_PACKET_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `R29/P36/N39`, the narrowest route-scoped local blocker on the attacked
source side was:

```text
explicit assignment witness of m2_psi1 to mu_m2_plus3_segment_psi1_psi4
```

`R30` does not pretend to prove that source-slot assignment witness.

It attacks only the next honest subobject:

```text
can that one source-slot assignment witness be split exactly into source
action-role and source eom-role assignment witnesses using the already
exported R26 role support?
```

This keeps the light issue explicit:

1. the shared kernel/light-facing channel remains already closed by `R14`,
2. `R30` touches only one non-light direct `m2` source slot on the already
   exported sufficient route.

## Inputs reused

1. `R29`
   - exact one-pair slotwise split of the combined assignment witness.
2. `R26`
   - exact source action role `m2_psi1*psi1**2/2`,
   - exact source eom role `m2_psi1*psi1(x)`.

## Result of `R30`

`R30` materializes one exact route-scoped role split:

1. the missing source-slot witness under attack is
   `explicit_assignment_witness_of_m2_psi1_to_mu_m2_plus3_segment_psi1_psi4`,
2. on this route, that source-slot witness can only appear through the two
   narrower role-specific witnesses:
   `explicit_source_action_role_assignment_witness_for_m2_psi1_to_mu_m2_plus3_segment_psi1_psi4`,
   `explicit_source_eom_role_assignment_witness_for_m2_psi1_to_mu_m2_plus3_segment_psi1_psi4`,
3. `R30` does **not** claim that either role-specific witness is present.

So, on this one source slot only, the repo now exports not only the slotwise
assignment gap but also the exact current reason why it still fails:
neither source-role assignment witness is present.

## Honest frontier after `R30`

On the direct formal coefficient-family route, the host route is reduced to:

1. explicit zero witness for direct quartic-like `g4` family `c1s1`
   shift defect,
2. explicit zero witness for direct quintic-like `g6` family `c1s1`
   shift defect,
3. explicit zero witness for direct yukawa-like `gY` family `c1s1`
   shift defect,
4. explicit source action-role assignment witness for `m2_psi1` to
   `mu_m2_plus3_segment_psi1_psi4`,
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

## What `R30` does not claim

`R30` does not claim:

- theorem-level PASS,
- full-closure PASS,
- that `m2_psi1 = m2_psi4`,
- that any common plus3 carrier-segment parameter actually exists,
- that either source-role assignment witness is present,
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

1. rerun the direct formal family route after this exact source-role split,
2. accept only:
   - a shorter direct-route missing-object list,
   - or the unchanged negative route.
