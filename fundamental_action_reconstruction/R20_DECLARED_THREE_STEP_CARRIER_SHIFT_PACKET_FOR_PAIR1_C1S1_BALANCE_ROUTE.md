# R20 Declared Three-Step Carrier Shift Packet For Pair1 C1S1 Balance Route

Status: `R20_EXECUTED_DECLARED_THREE_STEP_CARRIER_SHIFT_PACKET_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `R19/P26/N29`, the narrowest remaining local blocker on the host route
was:

```text
explicit balance witness for the declared pair1 c1s1 equation
```

`R20` does not pretend to prove that balance witness.

It asks the narrower constructive question:

```text
does the repo already contain an exact declared carrier shift that flips the
pair1 c1s1 sign and maps the positive support classes onto the negative ones?
```

This keeps the light issue explicit:

1. the shared kernel/light-facing channel remains already closed by `R14`,
2. `R20` touches only the non-light residual local diagonal complement.

## Inputs reused

1. `R11`
   - explicit declared control transport matrix on `psi0..psi11`.
2. `R18`
   - exact coefficient-class reduction of the declared `pair1` residual block.
3. `R19`
   - exact balance reduction of the declared `pair1` `c1s1` branch.

## Result of `R20`

`R20` materializes:

1. the exact declared carrier shift
   `S_{+3} : psi_i -> psi_{i+3 mod 12}`,
2. the exact transported pair1 action
   `S_{+3}(c1) = s1`,
   `S_{+3}(s1) = -c1`,
3. the exact support-class map
   `Sigma_psi1_psi7 -> Sigma_psi4_psi10`,
   `Sigma_psi2_psi8 -> Sigma_psi5_psi11`,
4. the equivalent reduction:
   the missing `c1s1` balance witness is reduced to one declared
   `S_{+3}`-equivariance witness on the positive support sum.

So the old missing object

```text
explicit balance witness for the declared pair1 c1s1 equation
```

is now reduced to the narrower object

```text
explicit declared plus3 shift-equivariance witness for the pair1 c1s1 support
sum
```

## Honest frontier after `R20`

After `R20`, the host route is reduced to:

1. explicit declared `S_{+3}` shift-equivariance witness for the `pair1`
   `c1s1` support sum,
2. explicit zero witness for the declared `pair1` `c1c1` equation,
3. explicit zero witness for the declared `pair1` `s1s1` equation,
4. `QW-2191` physical canonicalization boundary.

The already closed light-facing kernel channel remains unchanged.

## What `R20` does not claim

`R20` does not claim:

- theorem-level PASS,
- full-closure PASS,
- that the declared `S_{+3}` shift is physically canonical,
- that the declared `pair1` `c1s1` support sum is shift-invariant,
- that the declared `pair1` `c1c1` equation holds,
- that the declared `pair1` `s1s1` equation holds,
- that `QW-2191` is discharged,
- that selector closure is obtained,
- that ToE is closed.

## Recommended next move

The correct next move is now:

1. rerun the host-matching route after this declared shift packet,
2. accept only:
   - a shorter missing-object list,
   - or the unchanged negative route.
