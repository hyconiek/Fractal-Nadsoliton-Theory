# P2023 S973 Strict Cutkosky BRST Projector Source-Status Theorem

Status: `P2023_BRST_PROJECTOR_SOURCE_STATUS_LOCAL_TRANSVERSE_ONLY_WITH_TRACE_NO_FALSE_PASS`
As of: `2026-05-19`

## Goal

P2022 identified the remaining unitarity blocker as BRST/ghost/scheme data, not
the absence of the minimal tree `hAA` source chain.  P2023 performs the next
strict check: separate the already available local transverse polarization
projector from a genuine channel-level BRST physical-state projector.

## Result

P2023 credits the existing local result:

```text
P1956/P2020 provide the local on-shell transverse gauge-gauge projector and the
matching P2020 real {plus,cross} CutSum matrix.
```

But it does not promote that result to BRST cohomology.  Four theorem-level
objects remain missing:

1. a nilpotent `Q_BRST` operator map on the same channel Hilbert space,
2. a BRST physical-state/cohomology projector for `graviton -> gauge_gauge`,
3. a ghost-sector exclusion/cancellation trace in the same phase-space
   normalization,
4. a same-scheme bridge from those BRST data to `DiscM_common_basis`.

## Symbolic gap witness

The P2020 no-identical-symmetry local transverse trace is

```text
Trace(CutSum_transverse) = 2/pi.
```

P2023 records a formal physical trace

```text
Trace_BRST_phys = 2/pi - Q_P_T - GhostTrace + BRST_exact.
```

Without exported values for `Q_P_T`, `GhostTrace`, and `BRST_exact`, the optical
defect is underdetermined.  The optimistic assignment gives defect `0`, while
nonzero ghost or BRST-exact corrections give nonzero defects.  Therefore the
local transverse projector cannot be read as a BRST-projected optical theorem.

## No-false-pass boundary

P2023 is not:

1. a BRST nilpotency theorem,
2. a BRST cohomology projector export,
3. a ghost cancellation proof,
4. a `DiscM = CutSum` theorem,
5. all-state unitarity,
6. `QW-2191` discharge,
7. ToE closure.

It is a source-status theorem that protects the valid P1956/P2020 local
projector from being over-interpreted.

## Progress toward ToE

This is progress by sharpening the unitarity frontier.  The project no longer
confuses two different objects:

```text
local transverse polarization projector != BRST physical-state projector.
```

The exact tree CutSum remains useful, but theorem-grade unitarity now requires
an explicit BRST/ghost layer in the same scheme.

## Next honest step

Build P2024 by exporting the smallest channel-level BRST data object: `Q_BRST`
action on gauge, ghost, and antighost cut states; a machine-checkable
nilpotency condition `Q_BRST^2=0`; and a ghost-sector trace convention compatible
with the P2020 phase-space normalization.
