# R70 Direct M2 Psi5 Target EOM Common Monomial Support Packet

Status: `R70_EXECUTED_DIRECT_M2_PSI5_TARGET_EOM_COMMON_MONOMIAL_SUPPORT_PACKET_NO_FALSE_PASS`  
As of: `2026-03-16`

## Goal

After `R63`, the narrowest route-scoped blocker on the attacked target eom side
is:

```text
explicit target eom-role assignment witness for m2_psi5 to
mu_m2_plus3_segment_psi2_psi5
```

`R70` does not pretend to prove that eom-role assignment witness.

It attacks only the next honest subobject:

```text
export the fixed common local support psi5(x) and reduce the eom-role
assignment witness to one coefficient-identification witness on that support
```

## Inputs reused

1. `R63`
   - exact target eom term `m2_psi5*psi5(x)`,
   - declared lifted target eom term `mu_m2_plus3_segment_psi2_psi5*psi5(x)`.

## Result of `R70`

`R70` fixes the common target-eom support `psi5(x)` and reduces the missing
eom-role assignment witness to the narrower coefficient-identification gap:

```text
m2_psi5 vs mu_m2_plus3_segment_psi2_psi5 on common psi5(x) support
```

## What `R70` does not claim

`R70` does not claim:

- theorem-level PASS,
- full-closure PASS,
- that `m2_psi5 = mu_m2_plus3_segment_psi2_psi5`,
- that `m2_psi2 = m2_psi5`,
- that `psi5(x)` may be divided out or treated as a nonzero factor,
- that any direct `g4/g6/gY` family defect vanishes,
- that any `pair1` residual zero equation holds,
- that `QW-2191` is discharged,
- that selector closure is obtained,
- that ToE is closed.

## Recommended next move

The correct next move is now:

1. materialize the exact coefficient defect polynomial on this common support,
2. keep the missing object explicitly as a zero witness (no false PASS).

