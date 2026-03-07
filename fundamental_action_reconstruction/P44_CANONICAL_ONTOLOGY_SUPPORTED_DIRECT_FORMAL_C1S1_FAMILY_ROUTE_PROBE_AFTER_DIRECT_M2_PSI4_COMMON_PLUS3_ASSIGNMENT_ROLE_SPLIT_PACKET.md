# P44 Canonical Ontology Supported Direct Formal C1S1 Family Route Probe After Direct M2 Psi4 Common Plus3 Assignment Role Split Packet

Status: `P44_EXECUTED_CANONICAL_ONTOLOGY_SUPPORTED_DIRECT_FORMAL_C1S1_FAMILY_ROUTE_AFTER_R35_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `AX10/AX11/P43/N46`, the attacked target-side blocker is:

```text
explicit assignment witness of m2_psi4 to mu_m2_plus3_segment_psi1_psi4
```

`R35` now adds one honest local reduction packet:

```text
exact target action/eom role split of that one target-slot assignment witness
```

`P44` reruns the canonical-ontology-supported direct family route after that
addition.

## Result

The route remains negative:

```text
CANONICAL_ONTOLOGY_SUPPORTED_SOURCE_SIDE_CLOSED_BUT_TARGET_M2_PSI4_ROLE_ASSIGNMENT_GAPS_REMAIN_AFTER_R35
```

## Why it still fails

`R35` gives only:

1. exact route-scoped role splitting of the one target-slot assignment witness,
2. two explicit target-role missing witnesses:
   `m2_psi4*psi4**2/2 -> mu_m2_plus3_segment_psi1_psi4*psi4**2/2`,
   `m2_psi4*psi4(x) -> mu_m2_plus3_segment_psi1_psi4*psi4(x)`.

It still does **not** give:

1. either actual target-role assignment witness,
2. the other three direct `m2` pairwise witnesses,
3. the direct `g4/g6/gY` zero witnesses,
4. the declared `pair1` `c1c1` zero equation,
5. the declared `pair1` `s1s1` zero equation,
6. selector-relevant canonicalization beyond `QW-2191`.

## Real reduction after `R35`

`P44` does not claim that the main `R21/P28` host frontier is globally solved.

It claims only this narrower route-scoped decomposition:

```text
explicit assignment witness of m2_psi4 to mu_m2_plus3_segment_psi1_psi4
  -> explicit target action-role assignment witness
  + explicit target eom-role assignment witness
```

So the route is shorter by one local layer on the attacked target slot, but
still negative.

## What `P44` does not claim

`P44` does not claim:

- theorem-level PASS,
- full-closure PASS,
- that `m2_psi1 = m2_psi4`,
- that any common plus3 carrier-segment parameter actually exists,
- that either target-role assignment witness is present,
- that any other direct `m2` pairwise equality holds,
- that the direct `m2` shift-equivariance holds,
- that any direct `g4/g6/gY` family defect vanishes,
- that the declared `pair1` `c1c1` equation holds,
- that the declared `pair1` `s1s1` equation holds,
- that `QW-2191` is discharged,
- that the direct route closes ToE.

## Recommended next move

The correct next move is now:

1. either attack one of the two target-role assignment witnesses for
   `m2_psi4`,
2. or separately attack one of the remaining direct `m2` pairwise witnesses,
3. or separately attack one of the direct `g4/g6/gY` zero witnesses,
4. or keep both the main host route and this direct formal family route
   negative.
