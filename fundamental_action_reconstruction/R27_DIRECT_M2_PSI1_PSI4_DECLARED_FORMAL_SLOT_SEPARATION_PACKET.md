# R27 Direct M2 Psi1 Psi4 Declared Formal Slot Separation Packet

Status: `R27_EXECUTED_DIRECT_M2_PSI1_PSI4_DECLARED_FORMAL_SLOT_SEPARATION_PACKET_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `R26/P33/N36`, the narrowest route-scoped local blocker on the single
direct `m2` pair was:

```text
explicit coefficient-identification witness for the declared role-matched
direct m2 pair m2_psi1 = m2_psi4
```

`R27` does not pretend to prove that witness.

It attacks only the next honest subobject:

```text
can the single coefficient-identification gap be reduced to one narrower
common-parameter-source or symbol-identification witness between two explicit
named slots of the same exported m2 family?
```

This keeps the light issue explicit:

1. the shared kernel/light-facing channel remains already closed by `R14`,
2. `R27` touches only one non-light direct `m2` pair on the already exported
   sufficient route.

## Inputs reused

1. `R26`
   - exact one-pair action/eom role match for `m2_psi1 / m2_psi4`.
2. `QW-2164`
   - canonical continuum potential with the full `m2_psi_i` family.
3. `QW-2165`
   - full canonical lagrangian density and local eom exports.

## Result of `R27`

`R27` materializes one route-scoped declared formal slot-separation packet:

1. current canonical export carries the family
   `m2_psi0, ..., m2_psi11`,
2. `m2_psi1` and `m2_psi4` are two distinct named slots in that same family,
3. `R26` already gives exact role matching for that pair under the declared
   `+3` shift,
4. therefore the remaining one-pair gap is narrowed to:
   `explicit_common_parameter_source_or_symbol_identification_witness_for_the_declared_formal_m2_slots_m2_psi1_and_m2_psi4`.

So, on this one pair only, the repo now exports not only role matching but also
the exact current export-level reason why equality still does not follow:
the two coefficients still live in two distinct named formal slots.

## Honest frontier after `R27`

On the direct formal coefficient-family route, the host route is reduced to:

1. explicit zero witness for direct quartic-like `g4` family `c1s1`
   shift defect,
2. explicit zero witness for direct quintic-like `g6` family `c1s1`
   shift defect,
3. explicit zero witness for direct yukawa-like `gY` family `c1s1`
   shift defect,
4. explicit common-parameter-source or symbol-identification witness for the
   declared formal direct `m2` slots `m2_psi1` and `m2_psi4`,
5. explicit pairwise matching witness for `m2_psi7 = m2_psi10`,
6. explicit pairwise matching witness for `m2_psi2 = m2_psi5`,
7. explicit pairwise matching witness for `m2_psi8 = m2_psi11`,
8. explicit zero witness for the declared `pair1` `c1c1` equation,
9. explicit zero witness for the declared `pair1` `s1s1` equation,
10. `QW-2191` physical canonicalization boundary.

The already closed light-facing kernel channel remains unchanged.

## What `R27` does not claim

`R27` does not claim:

- theorem-level PASS,
- full-closure PASS,
- that `m2_psi1 = m2_psi4`,
- that the two named slots are already physically identified,
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

1. rerun the direct formal family route after this one-pair formal
   slot-separation packet,
2. accept only:
   - a shorter direct-route missing-object list,
   - or the unchanged negative route.
