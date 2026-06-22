# P3029/S1979 matter spectral observable generator obstruction

Status: `P3029_MATTER_SPECTRAL_OBSERVABLE_GENERATOR_OBSTRUCTION_NO_CLASSICAL_EXPORT`

## Finite certificate
- signature length: `12`
- U(12) rows: `4`
- invariant rows: `4`
- observer-independent generator accepted: `True`
- matter-sector export accepted: `False`

## Decision
A real observer-independent observable generator for the P3028 matter_fields row was constructed: the sorted DFT magnitude signature of K_strict_gate on Z12.  It has explicit types and is invariant under U(12) relabeling, but it is not a matter-sector export because it lacks a field representation/localizer, mass/coupling provenance, and selector/sector source.

## Recommendation
Do not promote spectral magnitude signatures to matter physics.  The next proof-grade move may attack exactly one missing matter atom for this generator: a field-representation/localizer theorem or a mass/coupling provenance theorem; otherwise return to the P3028 lattice and choose another single foundation atom.
