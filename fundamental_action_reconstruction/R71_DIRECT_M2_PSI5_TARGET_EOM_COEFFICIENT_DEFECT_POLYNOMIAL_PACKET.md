# R71 Direct M2 Psi5 Target EOM Coefficient Defect Polynomial Packet

Status: `R71_EXECUTED_DIRECT_M2_PSI5_TARGET_EOM_COEFFICIENT_DEFECT_POLYNOMIAL_PACKET_NO_FALSE_PASS`  
As of: `2026-03-16`

## Goal

After `R70`, the narrowest route-scoped blocker on the attacked target eom side
is:

```text
explicit target eom monomial coefficient-identification witness for
m2_psi5 and mu_m2_plus3_segment_psi2_psi5 on common psi5(x) support
```

`R71` does not pretend to prove that coefficient-identification witness.

It attacks only the next honest subobject:

```text
materialize the exact coefficient defect polynomial whose vanishing would
imply that missing witness, without dividing by psi5(x) and without any
nonzero-factor claim
```

## Inputs reused

1. `R70`
   - exact target eom term on `psi5(x)`,
   - declared lifted target eom term on the same support.

## Result of `R71`

`R71` exports the exact target-eom coefficient defect polynomial:

```text
(m2_psi5) - (mu_m2_plus3_segment_psi2_psi5)
```

and records the exact local eom defect expression on common support:

```text
((m2_psi5) - (mu_m2_plus3_segment_psi2_psi5))*psi5(x)
```

The remaining missing object is the corresponding explicit zero witness.

## What `R71` does not claim

`R71` does not claim:

- theorem-level PASS,
- full-closure PASS,
- that the coefficient defect polynomial vanishes,
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

1. integrate this packet into the tracked direct-formal route probe,
2. keep the missing object explicitly as a zero witness (no false PASS).

