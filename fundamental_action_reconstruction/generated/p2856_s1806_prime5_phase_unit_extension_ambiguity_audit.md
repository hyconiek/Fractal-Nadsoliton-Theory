# P2856/S1806 prime-5 phase-unit extension ambiguity audit

Status: `P2856_PRIME5_PHASE_UNIT_EXTENSION_AMBIGUITY_AUDIT_NO_CLOSURE`

## Exact representation after importing prime 5
- omega=743/4000; first exact denominators=[4000]
- phi=13/80; first exact denominators=[80, 160, 240, 320, 400, 480, 640, 720, 800, 960, 1200, 1280]

## Local ambiguity witnesses
- local window: `1/200` around the strict tuple
- witness count reported: `12`
- omega=3/16, phi=13/80, L1 error=7/4000
- omega=3/16, phi=329/2025, L1 error=577/324000
- omega=3/16, phi=79/486, L1 error=1751/972000
- omega=3/16, phi=508/3125, L1 error=181/100000
- omega=3/16, phi=499/3072, L1 error=697/384000

## Boundary
P2856 grants the missing prime-5 denominator support and confirms exact representation, but finite local enumeration finds multiple non-strict extended-lattice pairs with the same Z12 phase-bit profile.  The prime-5 extension therefore supplies coordinate capacity, not a non-premise source law selecting omega=743/4000 and phi=13/80.

## Recommendation
Do not replay pure Z12 lattices, prime-5 representability, or phase-bit equivalence as a source.  The next proof-grade move requires a genuinely new source-selection law for the prime-5 phase unit and the exact omega/phi numerators, or a genuinely new eta/beta source law.  Without that new typed premise, preserve the no-new-live-frontier certificate.
