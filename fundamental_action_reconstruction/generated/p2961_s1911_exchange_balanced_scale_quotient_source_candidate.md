# P2961/S1911 exchange-balanced scale quotient source candidate

Status: `P2961_EXCHANGE_BALANCED_SCALE_QUOTIENT_SOURCE_CANDIDATE_NO_STRICT_EXPORT`

## Quotient certificate
- bounded pair rows: `144`
- exchange orbit count: `46`
- fixed orbits: `[{'orbit': [[1, 1], [1, 1]], 'representative': [1, 1], 'is_fixed_orbit': True, 'selected_by_exchange_balance': True, 'vector': [1, 2, 2, 2, 2], 'sum': 9, 'equals_target': True}]`
- unique fixed orbit is target: `True`
- target sum from fixed orbit: `9`
- developmental source candidate exported: `True`
- strict nadsoliton exchange-symmetry source exported: `False`
- unit-bearing nonproxy coupling exported: `False`
- strict ratio-package source exported: `False`
- acceptance matrix rows/developmental/strict: `64/4/1`

## Lay summary
P2961 upgrades one P2960 candidate: the scale quotient and exchange-fixed primitive orbit are now explicit, and the unique fixed orbit produces the exact P2938 vector and sum 9 without using a target_sum=9 cut.
This is a real source-candidate theorem, not full strict closure.  The missing step is to prove that the nadsoliton itself exports the K/C exchange symmetry and then couple the selected quotient to a unit-bearing nonproxy action density.

## Recommendation
Do not replay exchange-fixed quotient enumeration, target_sum=9 cuts, primitive equal-summand convention, bounded localizer variants, K+C decompositions, beta-scale normalization, scalar Euler insertion, or count/role aliases.  The next proof-grade move must either prove a strict nadsoliton source theorem for the K/C exchange symmetry used by P2961, or construct the unit-bearing nonproxy action-density coupling receiving the P2961 quotient-selected vector; otherwise pivot outside the ratio-package lane while preserving the P2929-P2961 no-strict-export boundary.
