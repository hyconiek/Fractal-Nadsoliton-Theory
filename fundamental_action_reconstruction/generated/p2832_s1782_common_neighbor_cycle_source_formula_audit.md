# P2832/S1782 common-neighbor/cycle source formula audit

Status: `P2832_COMMON_NEIGHBOR_CYCLE_SOURCE_FORMULA_NO_GO_NO_COUPLING_NO_CLOSURE`

## Candidate formula
Q_cyc(G)=(triangles(G), 4-cycles(G), histogram_{pairs {u,v}} |N(u)∩N(v)|), normalized by fixed n=16 and degree=4 constants when used as a dimensionless coefficient vector

## Finite counts
- decoded_graph_count=16828
- candidate_class_count=344
- candidate_collision_class_count=272
- candidate_collision_graph_count=16756
- candidate_max_class_size=788
- candidate_defect_after_formula=16484

## Acceptance
- accepted_as_source_law_coupling=False
- accepted_as_bounded_formula_no_go=True

## Boundary
P2832 tests one genuinely non-label graph formula after P2831: the exact common-neighbor/cycle profile Q_cyc.  The formula is concrete and dimensionlessly normalizable on the fixed 16-node 4-regular carrier, but the finite computation yields only 344 classes, with 272 residual collision classes covering 16,756 graphs and max class size 788.  Since it does not separate the P2830 carrier and still exports no variational derivative or typed graph-to-K/L_total coupling theorem, it is rejected as a source-law/coupling promotion.

## Recommendation
Do not promote low-order common-neighbor/cycle profiles to dynamics.  The next admissible source-law move must add one higher-order non-label formula that is not merely the digest label, preferably an exact edge-toggle response polynomial or typed action-density expression with a proved variational derivative into K or L_total; otherwise pivot away from graph-source coupling and preserve the P2831-P2832 no-coupling boundary.
