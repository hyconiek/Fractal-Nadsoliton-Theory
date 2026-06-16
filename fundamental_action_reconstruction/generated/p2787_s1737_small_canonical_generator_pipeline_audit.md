# P2787/S1737 small canonical generator pipeline audit

Status: `P2787_SMALL_CANONICAL_GENERATOR_PIPELINE_AUDIT_NO_CLOSURE`

## Exact small-pipeline result
- labeled_candidate_count=19355
- connected_labeled_candidate_count=19355
- isomorphism_class_count=6
- pair_count_after_quotient=15
- exact_charpoly_collision_counts={'adjacency_charpoly_coefficients': 0, 'laplacian_charpoly_coefficients': 0, 'signless_laplacian_charpoly_coefficients': 0}
- all_pairs_separated_by_all_three_exact_charpolys=True

## Decision
The exact generator/quotient/charpoly pipeline is validated on the complete 8-node 4-regular class, but it is not the required full 16-node class and exports no strict K/L_total source law.

## Recommendation
Use P2787 only as a complete small-class pipeline validation.  The next honest move is exactly one of: supply/import an actual certified full connected 16-node 4-regular generator artifact/toolchain with graph6/hash provenance and run the same exact quotient/charpoly audit there; or export a strict nadsoliton spectral action/source law fixing the admissible class, target spectrum, and K/L_total coupling before testing.  Otherwise preserve the P2697-P2787 no-canonical-geometry/no-closure certificate.
