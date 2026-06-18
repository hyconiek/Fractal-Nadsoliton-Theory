# P2867/S1817 coupled Z12-smooth weighted log-functional no-source audit

Status: `P2867_COUPLED_Z12_SMOOTH_WEIGHTED_LOG_FUNCTIONAL_NO_SOURCE_AUDIT_NO_CLOSURE`

## Coupled weighted-functional scan
- candidate class: `F(w)=sum_{d=1}^{11} w_d log(d), with every w_d in Q whose denominator support is contained in {2,3}`
- coordinate forcing: `{'prime11_contributors': [11], 'forced_w_11': {'fraction': '9/5', 'numerator': 9, 'denominator': 5, 'denominator_prime_factors': {5: 1}, 'z12_compatible_denominator': False}, 'forced_w_11_denominator_primes': [5], 'forced_w_11_is_z12_smooth': False, 'obstruction': 'Since only d=11 contributes to the prime-11 coordinate, any exact weighted log functional must set w_11=9/5, which imports denominator prime 5 outside the Z12-smooth denominator support {2,3}.'}`
- bounded scan: `{'denominators_scanned': [1, 2, 3, 4, 6, 8, 9, 12], 'numerator_abs_max': 24, 'candidate_count': 217, 'exact_w11_hits': [], 'best_w11': {'fraction': '16/9', 'numerator': 16, 'denominator': 9, 'denominator_prime_factors': {3: 2}, 'z12_compatible_denominator': True}, 'best_prime11_coordinate_error': {'fraction': '-1/45', 'numerator': -1, 'denominator': 45, 'denominator_prime_factors': {3: 2, 5: 1}, 'z12_compatible_denominator': False}, 'exact_w11_absent': True}`

## Boundary
P2867 tests a genuinely coupled weighted-log functional rather than pairwise recombination.  Exact equality forces w_11=9/5 by the prime-11 coordinate, but that forced weight imports denominator prime 5 and is outside the pure Z12-smooth coefficient/source class.  The imported exact witness remains representability, not sourcehood.

## Recommendation
Do not replay Z12-smooth weighted log functionals or imported w_11=9/5 as boundary sourcehood.  A next proof-grade move must introduce a new strict mechanism that sources denominator prime 5 and the d=11 endpoint together, or pivot to a different new typed object; otherwise preserve no-new-live-frontier.
