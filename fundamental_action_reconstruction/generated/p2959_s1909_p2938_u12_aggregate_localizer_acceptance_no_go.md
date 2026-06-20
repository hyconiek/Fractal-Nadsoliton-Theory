# P2959/S1909 P2938 U(12) aggregate localizer acceptance no-go

Status: `P2959_P2938_U12_AGGREGATE_LOCALIZER_ACCEPTANCE_NO_GO_NO_STRICT_EXPORT`

## Localizer certificate
- candidate lattice rows: `15`
- exact target pair: `[1, 1]`
- exact target vector: `[1, 2, 2, 2, 2]`
- predicates selecting exact target only: `['primitive_equal_summand', 'positive_and_primitive_equal_summand']`
- all exact-selecting predicates strict-exported: `False`
- predicate not target-coded or conventional: `False`
- strict nadsoliton functor exports predicate: `False`
- strict equal-weight source theorem exported: `False`
- downstream beta/unit coupling exported: `False`
- P2958 functor/localizer obligation discharged: `False`
- acceptance matrix rows/accepted: `64/1`

## Boundary
P2959 attacks the missing P2958 nadsoliton-to-U12 functor/localizer.  In the bounded localizer lattice, sum=9 is target-coded and nonunique, while the target pair (a,b)=(1,1) is selectable by primitive equal-summand predicates; none is strict unless a nadsoliton functor exports it.  Current artifacts export no such functor, no strict equal-weight theorem, and no beta/unit coupling.

## Recommendation
Do not continue bounded localizer-predicate variants, K+C decompositions, symbolic weight-family variants, beta-scale normalization, P2601 replay, scalar Euler insertion, or count/role aliases.  A next proof-grade move must either introduce a genuinely new strict nadsoliton localizer law not equivalent to target-coded sum=9 cuts or primitive equal-weight convention, construct a unit-bearing nonproxy coupling outside this localizer lane, or pivot outside the ratio-package lane while preserving the P2929-P2959 no-strict-export boundary.
