# P2747/S1697 Z12 cubic phase orbit signed observable audit

Status: `P2747_CUBIC_PHASE_POLARITY_SOURCE_GAP`

## Finite cubic-phase audit
- coefficient_triple_count=1728
- pointwise signs: +396 / -396 / 0=936
- global_signed_sum=0
- affine_orbit_count=180
- affine_orbit_sizes=[1, 2, 3, 4, 6, 8, 12, 16, 24, 48]
- nonzero_orbit_coefficient_count=44
- positive_nonzero_coefficients=22
- negative_nonzero_coefficients=22
- nonzero_orbit_coefficient_histogram={'-8': 4, '-4': 9, '-2': 7, '-1': 2, '1': 2, '2': 7, '4': 9, '8': 4}

## Theorem statement
Cubic Z12 phase sums are a real pivot away from the P2746 Gauss-selector gap: 1728 coefficient triples produce 396 positive, 396 negative, and 936 zero imaginary signs.  Affine quotienting gives 180 coefficient orbits, with 44 nonzero signed-sum coefficients.  The nonzero coefficients are globally polarity-balanced by counts and values (-8/-4/-2/-1 paired with +1/+2/+4/+8), so the object supplies orbit-safe signed coefficients but still no strict source selecting one orbit/polarity or coupling it to P2721.

## Recommendation
Do not promote the P2747 cubic phase observable to lambda/P2721 fixing.  It gives many nonzero affine-orbit signed coefficients, but the coefficient family is still polarity-balanced and lacks a strict orbit/polarity source plus P2721 coupling theorem.  The next proof-grade move should audit exactly one missing premise: a strict cubic coefficient-orbit/polarity selector law with explicit P2721 coupling; if no such law is available, pivot to a different typed object or preserve the P2697-P2747 no-new-live-frontier certificate.
