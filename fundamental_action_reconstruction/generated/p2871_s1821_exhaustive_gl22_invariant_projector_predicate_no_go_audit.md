# P2871/S1821 exhaustive GL(2,2)-invariant projector-predicate no-go audit

Status: `P2871_EXHAUSTIVE_GL22_INVARIANT_PROJECTOR_PREDICATE_NO_GO_AUDIT_NO_CLOSURE`

## Exhaustive predicate audit
- candidate class: `all Boolean selector predicates on the four Aut(Z12)-character point projectors, filtered by V4/GL(2,2) invariance`
- predicate count: `16`
- invariant predicates: `[{'selected': [], 'cardinality': 0, 'is_gl22_invariant': True, 'orbit_size': 1, 'orbit': [[]], 'selects_target_singleton_11': False}, {'selected': [1], 'cardinality': 1, 'is_gl22_invariant': True, 'orbit_size': 1, 'orbit': [[1]], 'selects_target_singleton_11': False}, {'selected': [5, 7, 11], 'cardinality': 3, 'is_gl22_invariant': True, 'orbit_size': 1, 'orbit': [[11, 5, 7]], 'selects_target_singleton_11': False}, {'selected': [1, 5, 7, 11], 'cardinality': 4, 'is_gl22_invariant': True, 'orbit_size': 1, 'orbit': [[1, 11, 5, 7]], 'selects_target_singleton_11': False}]`
- target singleton record: `{'selected': [11], 'cardinality': 1, 'is_gl22_invariant': False, 'orbit_size': 3, 'orbit': [[5], [7], [11]], 'selects_target_singleton_11': True}`

## Boundary
P2871 exhausts all Boolean selector predicates on the four character point projectors.  GL(2,2)-invariant predicates cannot select singleton 11; the singleton target predicate is a non-invariant endpoint-label import.

## Recommendation
Do not replay Boolean projector predicates or target singleton predicates as sourcehood.  A next proof-grade move must provide a new strict symmetry-breaking law whose own provenance fixes the nonidentity singleton and supplies the unit-bearing coefficient/coupling theorem, or pivot to a different new typed object; otherwise preserve no-new-live-frontier.
