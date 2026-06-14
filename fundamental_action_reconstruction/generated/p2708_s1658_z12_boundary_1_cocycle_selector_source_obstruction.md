# P2708/S1658 Z12 boundary 1-cocycle selector-source obstruction

Status: `P2708_Z12_BOUNDARY_1_COCYCLE_OBSTRUCTION_NO_STRICT_SELECTOR`

## Computation
- rank_d0=11
- h1_dimension=1
- invariant_subspace_dimension_under_inversion=0

## Candidate table
- `oriented_boundary_1_cocycle_line_H1_Z12`: strict_provider_obligation_met=False. H1 supplies ±omega, not a canonical choice of +omega.
- `Aut_Z12_invariant_orientation`: strict_provider_obligation_met=False. Full Aut(Z12) includes inversion, which kills every nonzero invariant orientation.
- `premise_signed_orientation`: strict_provider_obligation_met=False. This is premise-based sign selection, not a non-premise strict-sourced provider.

## Decision
The Z12 boundary 1-cocycle is a real new typed object with H1 dimension 1, but it exports only an orientation line ±omega.  Inversion in Aut(Z12) sends omega to -omega, so the Aut-invariant nonzero orientation subspace is zero-dimensional.  Therefore this object can support a premise orientation but cannot discharge QW-2191 or export a non-premise strict selector on current artifacts.

## Next honest step
A further admissible move must provide a strict source for the missing sign of the boundary 1-cocycle, or a different genuinely new typed object/provider outside the closed lanes.  Without that, preserve the P2697-P2708 no-new-live-frontier/no-strict-selector certificate.
