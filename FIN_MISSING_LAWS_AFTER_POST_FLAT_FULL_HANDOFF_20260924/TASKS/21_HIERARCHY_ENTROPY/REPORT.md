# Hierarchy entropy after quotienting identical phase labels

Status: **BOUND_AND_CORRECTION**

## Result
The labeled-tree count T_M=M!/2^(M-1) overcounts physically indistinguishable
hierarchies when many leaves share the same phase label. The exact orbit count
should be obtained by Burnside. A rigorous lower bound is T_M/prod n_a!, and a
fixed-q count-vector encoding gives an exponential (not superextensive) upper
bound. Thus the M log M label entropy need not survive after quotienting a
finite phase alphabet.

## Key formulas
\[T_M=M!/2^{M-1}.\]
For group G=prod S_{n_a}: \[N_{orbits}\ge T_M/|G|.\]
For fixed q, hierarchical node count-vectors give \[\log N_{orbits}=O(M).\]

## Caveat
Earlier informal simple upper bound by a multinomial coefficient was not certified; use Burnside/count-vector bounds.

## Next question
Compute exact orbit counts or entropy density for large finite phase ensembles.
