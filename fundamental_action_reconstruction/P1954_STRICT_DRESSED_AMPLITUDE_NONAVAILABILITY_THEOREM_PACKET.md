# P1954 Strict Dressed Amplitude Nonavailability Theorem Packet

Status: `OPEN_OBSTRUCTION_WITH_TRACE__DRESSED_AMPLITUDE_UNDERDETERMINED_NO_FALSE_PASS`
As of: `2026-05-17`

## Goal

After `P1953`, decide whether the current strict exports are enough to derive

```text
DressedCutkoskyAmplitude_graviton_to_gauge_gauge_strict_B1_v1
```

without inventing missing physics.

## Result

The local verdict is:

```text
FORMAL_NONAVAILABILITY_THEOREM_EXPORTED__DRESSED_AMPLITUDE_UNDERDETERMINED
```

This is not a no-go theorem about the theory. It is a theorem about the current
repository state.

## Proof Trace

The current exports are insufficient because:

1. `P1907` exports a sector-level non-skeleton `L_total` registry, but no 4D
   graviton-gauge-gauge Feynman vertex tensor.
2. `P1866` exports a 1D symbolic proxy chain, not a 4D expansion around
   `g = eta + kappa*h`.
3. `P1852` exports BRST/Cutkosky seed contracts, not a channel physical-state
   projector.
4. `P1862` explicitly says dressed residues are seed-inherited and still need
   full propagator computation.
5. `P1913` leaves the `grmix` `DiscM/CutSum` row as open symbolic placeholders
   under `MSbar_candidate`, not locked to `MSbar_B1_seed`.

Therefore the `P1953` dressed-amplitude interface cannot be satisfied from
current strict exports.

## Minimal Missing Data

The next solver must export at least:

```text
V_hAA_tensor_strict_B1_v1
BRSTPhysicalProjector_gauge_gauge_strict_B1_v1
DressedGravitonPropagatorPoleResidueTable_strict_B1_v1
DiscM_CutSum_grmix_MSbar_B1_seed_common_basis_v1
SchemeLock_MSbar_B1_seed_for_P1913_grmix_v1
```

## Safe Consequence

The following local results remain valid:

1. `P1950` declared B1 counterterm cancellation,
2. `P1951` seed phase-space positivity,
3. `P1952` QW-2049 seed-local positivity rectangle.

But they must not be promoted to:

1. full dressed Cutkosky equality,
2. BRST-projected optical theorem closure,
3. ghost-free dressed propagator theorem,
4. global `UR_link_theorem`.

## Outputs

- `p1954_s904_strict_dressed_amplitude_nonavailability_theorem.py`
- `generated/p1954_s904_strict_dressed_amplitude_nonavailability_theorem.json`

## Next Honest Step

Build `P1955` with high reasoning: derive
`V_hAA_tensor_strict_B1_v1` from a 4D metric perturbation expansion of
`L_total`, or prove that the current `L_total` export lacks the metric-density
detail needed for that derivation.
