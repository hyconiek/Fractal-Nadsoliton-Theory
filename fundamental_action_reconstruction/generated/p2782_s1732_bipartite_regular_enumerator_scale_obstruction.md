# P2782/S1732 bipartite regular enumerator scale obstruction

Status: `P2782_BIPARTITE_REGULAR_ENUMERATOR_SCALE_OBSTRUCTION_NO_CLOSURE`

## Exact DP count
- row_mask_count=70
- dp_cache_states=75465
- labeled_bipartite_regular_matrix_count=116963796250
- row_column_relabeling_group_size=1625702400
- naive_enumerator_blocked_in_repo=True

## Decision
The exact fixed-bipartition count is too large for a naive in-repo canonical enumerator, and no canonical-generation theorem/tool certificate is supplied.  No spectral source law or K/L_total coupling is exported.

## Recommendation
Do not start a naive full 16-node 4-regular enumeration inside the repo.  The next honest move is exactly one of: supply a canonical-generation theorem/tool certificate (for example a checked nauty/geng-style graph6 import with reproducible hashes) and then run the full-spectrum collision audit; or export a strict nadsoliton spectral action/source law that fixes the admissible class and target before testing.  Otherwise preserve the P2697-P2782 no-canonical-geometry/no-closure certificate.
