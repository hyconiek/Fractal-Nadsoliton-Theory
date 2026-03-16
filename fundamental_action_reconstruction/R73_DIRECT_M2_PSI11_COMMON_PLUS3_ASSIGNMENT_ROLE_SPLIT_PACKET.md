# R73 Direct M2 Psi11 Common Plus3 Assignment Role Split Packet

Status: `R73_EXECUTED_DIRECT_M2_PSI11_COMMON_PLUS3_ASSIGNMENT_ROLE_SPLIT_PACKET_NO_FALSE_PASS`  
As of: `2026-03-16`

## Goal

After `R56`, the narrowest route-scoped local blocker on the direct `psi8 / psi11`
plus3 segment is:

```text
explicit assignment witness of m2_psi11 to mu_m2_plus3_segment_psi8_psi11
```

`R73` does not pretend to prove that target-slot assignment witness.

It attacks only the next honest subobject:

```text
can that one target-slot assignment witness be split exactly into target
action-role and target eom-role assignment witnesses using the already
exported R50 role support?
```

This keeps the light issue explicit: the shared kernel/light-facing channel
remains already closed by `R14`, while `R73` touches only one non-light direct
`m2` target slot on the already exported sufficient route.

## Inputs reused

1. `R56`
   - exact slotwise split of the combined common-plus3 assignment witness.
2. `R50`
   - exact target action role `m2_psi11*psi11**2/2`,
   - exact target eom role `m2_psi11*psi11(x)`.

## Result of `R73`

`R73` materializes one exact route-scoped role split:

1. the missing target-slot witness under attack is
   `explicit_assignment_witness_of_m2_psi11_to_mu_m2_plus3_segment_psi8_psi11`,
2. on this route, that target-slot witness can only appear through the two
   narrower role-specific witnesses:
   `explicit_target_action_role_assignment_witness_for_m2_psi11_to_mu_m2_plus3_segment_psi8_psi11`,
   `explicit_target_eom_role_assignment_witness_for_m2_psi11_to_mu_m2_plus3_segment_psi8_psi11`,
3. `R73` does **not** claim that either role-specific witness is present.

## What `R73` does not claim

`R73` does not claim:

- theorem-level PASS,
- full-closure PASS,
- that `m2_psi8 = m2_psi11`,
- that any common plus3 carrier-segment parameter actually exists,
- that either target-role assignment witness is present,
- that any direct `g4/g6/gY` family defect vanishes,
- that any `pair1` residual zero equation holds,
- that `QW-2191` is discharged,
- that selector closure is obtained,
- that ToE is closed.

## Recommended next move

The correct next move is now:

1. reduce one of the two role-specific target assignment witnesses into a
   coefficient-level defect expression on fixed support (action or eom),
2. accept only an explicit missing-object list (no false PASS).

