# R44 Direct M2 Psi7 Common Plus3 Assignment Role Split Packet

Status: `R44_EXECUTED_DIRECT_M2_PSI7_COMMON_PLUS3_ASSIGNMENT_ROLE_SPLIT_PACKET_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `R43/P55/N58`, the narrowest route-scoped local blocker on the attacked
source side of the remaining `psi7 / psi10` pair is:

```text
explicit assignment witness of m2_psi7 to mu_m2_plus3_segment_psi7_psi10
```

`R44` does not pretend to prove that source-slot assignment witness.

It attacks only the next honest subobject:

```text
can that one source-slot assignment witness be split exactly into source
action-role and source eom-role assignment witnesses using the already
exported R40 role support?
```

This keeps the light issue explicit:

1. the shared kernel/light-facing channel remains already closed by `R14`,
2. the attacked source-side and attacked `m2_psi4` target-side blockers remain
   locally closed only on the canonical-ontology-supported lane from
   `AX10/AX11/AX12/AX13`,
3. `R44` touches only one non-light direct `m2` source slot on the already
   exported sufficient route.

## Inputs reused

1. `R43`
   - exact one-pair slotwise split of the combined assignment witness.
2. `R40`
   - exact source action role `m2_psi7*psi7**2/2`,
   - exact source eom role `m2_psi7*psi7(x)`.

## Result of `R44`

`R44` materializes one exact route-scoped role split:

1. the missing source-slot witness under attack is
   `explicit_assignment_witness_of_m2_psi7_to_mu_m2_plus3_segment_psi7_psi10`,
2. on this route, that source-slot witness can only appear through the two
   narrower role-specific witnesses:
   `explicit_source_action_role_assignment_witness_for_m2_psi7_to_mu_m2_plus3_segment_psi7_psi10`,
   `explicit_source_eom_role_assignment_witness_for_m2_psi7_to_mu_m2_plus3_segment_psi7_psi10`,
3. `R44` does **not** claim that either role-specific witness is present.

So, on this one source slot only, the repo now exports not only the slotwise
assignment gap but also the exact current reason why it still fails:
neither source-role assignment witness is present.

## Honest frontier after `R44`

On the canonical-ontology-supported direct formal coefficient-family route, the
frontier is now reduced to:

1. explicit zero witness for direct quartic-like `g4` family `c1s1`
   shift defect,
2. explicit zero witness for direct quintic-like `g6` family `c1s1`
   shift defect,
3. explicit zero witness for direct yukawa-like `gY` family `c1s1`
   shift defect,
4. explicit source action-role assignment witness for `m2_psi7` to
   `mu_m2_plus3_segment_psi7_psi10`,
5. explicit source eom-role assignment witness for `m2_psi7` to
   `mu_m2_plus3_segment_psi7_psi10`,
6. explicit assignment witness of `m2_psi10` to
   `mu_m2_plus3_segment_psi7_psi10`,
7. explicit pairwise matching witness for `m2_psi2 = m2_psi5`,
8. explicit pairwise matching witness for `m2_psi8 = m2_psi11`,
9. explicit zero witness for the declared `pair1` `c1c1` equation,
10. explicit zero witness for the declared `pair1` `s1s1` equation,
11. `QW-2191` physical canonicalization boundary.

The already closed light-facing kernel channel remains unchanged.

## What `R44` does not claim

`R44` does not claim:

- theorem-level PASS,
- full-closure PASS,
- that `m2_psi7 = m2_psi10`,
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

1. rerun the canonical-ontology-supported direct formal family route after this
   exact source-role split,
2. accept only:
   - a shorter route-local missing-object list,
   - or the unchanged negative route.
