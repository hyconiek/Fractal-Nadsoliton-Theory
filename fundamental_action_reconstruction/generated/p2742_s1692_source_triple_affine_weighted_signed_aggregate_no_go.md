# P2742/S1692 source-triple affine-weighted signed aggregate no-go

Status: `P2742_SOURCE_TRIPLE_AFFINE_WEIGHTED_SIGNED_AGGREGATE_NO_GO`

## Finite aggregate audit
- ordered_distinct_triples=1320
- affine_group_size=48
- affine_ordered_orbit_count=34
- affine_orbit_sizes=[24, 48]
- orbits_with_nonzero_signed_sum_coefficient=0
- signed_sum_linear_map_rank=0
- signed_sum_linear_map_nullity_over_orbit_weights=34
- one_hot_weight_crosscheck_nonzero_count=0
- all_unit_11_pairing_witnesses_pass=True

## Theorem statement
Every affine ordered orbit of distinct Z12 triples is paired by the inversion unit 11 into opposite cyclic chiralities, so each orbit has signed-sum coefficient zero.  Therefore any affine-invariant orbit weighting, over any coefficient field of characteristic not 2, has total signed aggregate zero.  The P2740 chirality cannot be made into a nonzero orbit-safe strict signed value by affine orbit weights after P2741 blocks fixed localizers.

## Recommendation
Do not continue the P2740/P2741 source-triple lane by affine orbit weighting: P2742 proves every affine ordered orbit has zero chirality signed-sum coefficient, so arbitrary affine-invariant weights still give total signed aggregate zero.  The next proof-grade move must pivot to a genuinely different strict signed observable with a computable nonzero orbit-safe signed value and a P2721 coupling theorem, or else preserve the P2697-P2742 no-new-live-frontier certificate.
