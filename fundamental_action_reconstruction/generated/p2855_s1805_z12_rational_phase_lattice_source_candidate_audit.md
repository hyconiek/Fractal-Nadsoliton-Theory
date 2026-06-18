# P2855/S1805 Z12 rational phase-lattice source candidate audit

Status: `P2855_Z12_RATIONAL_PHASE_LATTICE_SOURCE_CANDIDATE_AUDIT_NO_CLOSURE`

## Exact denominator audit
- omega=743/4000; factors={2: 5, 5: 3}; z12-compatible=False
- phi=13/80; factors={2: 4, 5: 1}; z12-compatible=False

## Best Z12-compatible approximations
- omega best=107/576; error=1/72000
- phi best=79/486; error=1/19440
- common pair omega=361/1944; phi=79/486; exact_pair=False

## Candidate matrix
- pure_z12_denominator_lattice_exact_source: finite_witness_passes=False; exports_strict_source_law=False; blocked: reduced strict denominators contain prime 5, outside the prime support of Z12.
- bounded_z12_lattice_best_approximation: finite_witness_passes=True; exports_strict_source_law=False; approximation exists but exact omega/phi sourcing fails; approximation is not a source law.
- common_z12_lattice_pair_approximation: finite_witness_passes=True; exports_strict_source_law=False; phase_bits_match_p2853=True; common-denominator approximation is nonexact; phase-bit agreement, if present, is only a coarse witness.
- import_prime5_phase_unit_extension: finite_witness_passes=True; exports_strict_source_law=False; allowing prime 5 represents the strict tuple exactly, but the prime-5 unit is imported rather than sourced from Z12.

## Boundary
P2855 tests one concrete candidate source class: a pure Z12-compatible rational phase lattice.  The strict tuple is exact rational data, but its reduced denominators contain prime 5, outside the Z12 denominator prime support {2,3}.  Enumerated Z12-compatible approximations are nonexact.  Importing a prime-5 unit would represent the tuple, but would be a new unsourced premise rather than a strict source theorem.

## Recommendation
Do not replay pure Z12 rational-lattice sourcing for omega/phi.  A next proof-grade move must either supply a genuine source for the prime-5 phase unit, or pivot to a genuinely new eta/beta source law.  Without one of those new typed premises, preserve the P2854 no-new-live-frontier certificate.
