# P2024 S974 Strict Cutkosky Minimal Finite BRST Quartet Model Theorem

Status: `P2024_PARTIAL_FINITE_BRST_QUARTET_MODEL_V11_EXPORTED_NOT_CHANNEL_PROJECTOR_DISCM_GHOST_CUT_STILL_OPEN_NO_FALSE_PASS`
As of: `2026-05-19`

## Goal

P2023 showed that the local transverse projector in P1956/P2020 must not be
promoted to a full BRST physical-state projector. P2024 performs the next
honest strict-side move: export the smallest finite BRST quartet model (v11 wording with nonuniqueness, symbolic-lock, scheme-tag lock, digest-selfcheck, python-major lock, Q^2 spectral-lock, and Q^2 characteristic-polynomial lock hardening) that
is actually machine-checkable, while keeping the optical theorem obstruction
open.

## Exported minimal object

P2024 exports a finite one-gauge-polarization-layer BRST quartet model with basis

```text
{ gauge_T1, gauge_T2, gauge_L, ghost, antighost, B }.
```

The BRST action is

```text
Q(gauge_T1)  = 0,
Q(gauge_T2)  = 0,
Q(gauge_L)   = ghost,
Q(ghost)     = 0,
Q(antighost) = B,
Q(B)         = 0.
```

The `gauge_T1/gauge_T2` labels are deliberately not the P2020 graviton
`plus/cross` labels. This object is anchored to the local gauge-sector BRST
differential already exported by P1961 and to the P1956 transverse gauge
projector.

## Exact checks

The exported matrix satisfies

```text
Q^2 = 0,
rank(Q) = 2,
dim ker(Q) = 4,
dim im(Q) = 2,
dim H_Q = dim ker(Q) - dim im(Q) = 2.
```

Thus the finite model has two gauge-transverse cohomology representatives. P2024(v5) additionally exports a transverse-label nonuniqueness witness: swapping gauge_T1/gauge_T2 leaves the BRST transverse subblock invariant, so no unique plus/cross identification can be claimed from this finite algebra alone. This
is a finite algebraic count only; it is not an amplitude-normalized identification
with the P2020 `plus/cross` CutSum matrix. A numpy/scipy cross-check
independently returns the same rank, kernel dimension, cohomology dimension, and
zero Frobenius norm for `Q^2`.

## Ghost-sector trace convention

P2024 also exports the minimal quartet supertrace convention

```text
weights(gauge_T1,gauge_T2,gauge_L,ghost,antighost,B) = (1,1,1,-1,-1,1).
```

The unphysical quartet has zero supertrace:

```text
Str(gauge_L,ghost,antighost,B) = 1 - 1 - 1 + 1 = 0,
```

leaving two gauge-transverse representatives at this finite projector layer.
This trace convention is not by itself an amplitude-normalized P2020 CutSum integral.

## Remaining Cutkosky gap

P2024(v8) also exports a same-scheme tag lock (`STRICT_P2020_PHASESPACE_SCHEME_V1`) to prevent semantic drift between the finite BRST layer and the declared P2020 phase-space normalization context.


P2024(v6) also exports a deterministic theorem-core digest over basis/Q-action/formal-gap/nonuniqueness flags, plus a self-consistency recomputation check, so accidental semantic drift is machine-detectable in regression checks.


This is not yet `DiscM = CutSum`. The symbolic post-P2024 obstruction is

```text
DiscM_loop - (AmpNorm*CutSum_tree + GhostCut_scheme + SchemeLock + WardLift
              + CohomologyAmplitudeBridge).
```

The optimistic unproved assignment can make this expression vanish, but a
nonzero ghost-cut or cohomology-to-amplitude bridge correction changes the
defect. Therefore the minimal quartet model narrows the obstruction; it does not
close unitarity.

## No-false-pass boundary

P2024 is not:

1. a full Hilbert-space BRST charge theorem,
2. an interacting BRST cohomology projector,
3. an identification of `gauge_T1/gauge_T2` with P2020 graviton `plus/cross`,
4. a ghost-cut amplitude derivation,
5. a cohomology-to-amplitude normalization bridge,
6. a same-scheme loop `DiscM_common_basis` computation,
7. a `DiscM = CutSum` theorem,
8. all-state unitarity,
9. `QW-2191` discharge,
10. ToE closure.

## Progress toward ToE

P2024 upgrades the unitarity route from “local transverse projector only” to an
executable finite BRST-quartet model. This is genuine progress because the route
now has a nilpotent algebraic BRST layer and a ghost-quartet trace convention.
The true bottleneck is now sharper: the missing same-scheme ghost-cut, Ward-lift,
cohomology-to-amplitude bridge, and loop-discontinuity data.

## Next honest step

Build P2025 by deriving the ghost-sector Cutkosky integrand, Ward-lift bridge,
and cohomology-to-amplitude normalization map in the same P2020 phase-space
normalization. Only after those data are exported should the project test a
same-scheme `DiscM_common_basis`.
