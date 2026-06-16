# P2798/S1748 external girth>=4 subtarget local girth-filter certificate

Status: `P2798_EXTERNAL_GIRTH4_SUBTARGET_LOCAL_GIRTH_FILTER_CERTIFICATE_NO_CLOSURE`

## External girth>=4 subtarget and local filter
- external_connected_16node_4regular_girth4_class_count=16828
- local_representative_count=8
- local_girth_at_least_4_count=6
- local_girth_at_least_4_labels=['circulant_pm1_pm3', 'circulant_pm1_pm4', 'circulant_pm1_pm6', 'circulant_pm1_pm7', 'torus_4x4', 'two_c8_layers_cross_pm0_pm4']
- local_triangle_containing_labels=['circulant_pm1_pm2', 'p2790_eighth_witness']
- external_girth4_gap_after_current_witnesses=16822
- local_girth4_coverage_fraction_exact=3/8414

## External source
- source=Markus Meringer Regular Graphs page
- url=https://www.mathe2.uni-bayreuth.de/markus/reggraphs.html
- linked_detail_url=https://www.mathe2.uni-bayreuth.de/markus/REGGRAPHS/16_4_4.html
- retrieved_utc_date=2026-06-16

## Decision
P2798 computes exact local girths and sizes a reachable external girth>=4 subtarget gap, but it does not import the shortcode graph list and therefore remains a target/gap certificate rather than generator provenance or a strict source law.

## Recommendation
Use P2798 only as a local exact-girth filter plus external girth>=4 target-gap certificate.  The next proof-grade move is to import the actual linked shortcode/graph-list artifact for the 16,828 girth>=4 graphs with hash/provenance and run exact quotient/charpoly/complement/orbit auditing, or else import the full 8,037,418-class toolchain/list, or export a strict spectral source law fixing the admissible class and K/L_total coupling.  Otherwise preserve the P2697-P2798 no-canonical-geometry/no-closure certificate.
